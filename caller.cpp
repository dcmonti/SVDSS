#include "caller.hpp"

#include <cctype>
#include <cstdlib>
#include <string>

void Caller::run() {
  config = Configuration::getInstance();

  // load reference genome and SFS
  load_chromosomes(config->reference);

  SFSData sfs_data = parse_sfsfile(config->sfs);
  SFSs = sfs_data.sfss;
  // SFS from FASTA have forward-strand coordinates that need inversion for
  // reverse reads. SFS from BAM are already in the right orientation.
  // UNKNOWN (legacy files without header) are treated as FASTA for backward
  // compatibility.
  bool sfs_from_fasta = (sfs_data.source != SFSSource::BAM);

  Clusterer C = Clusterer(&SFSs, sfs_from_fasta);
  C.run();

  spdlog::info("Calling SVs from {} clusters..", C.clusters.size());

  if (!config->normal_contigs_bam.empty()) {
    spdlog::info("Germline filter enabled: loading normal contigs from {}",
                 config->normal_contigs_bam);
    _p_normal_bam = new samFile *[config->threads];
    _p_normal_idx = new hts_idx_t *[config->threads];
    _p_normal_hdr = new bam_hdr_t *[config->threads];
    for (int i = 0; i < config->threads; i++) {
      _p_normal_bam[i] =
          hts_open(config->normal_contigs_bam.c_str(), "r");
      _p_normal_idx[i] = sam_index_load(_p_normal_bam[i],
                                         config->normal_contigs_bam.c_str());
      _p_normal_hdr[i] = sam_hdr_read(_p_normal_bam[i]);
    }
  }

  _p_svs.resize(config->threads);
  _p_alignments.resize(config->threads);
  _p_germline_regions.resize(config->threads);
  pcall(C.clusters);

  // NOTE: the normal-contigs BAM is kept open here; it is also used by the
  // germline filter on the clipped-SFS calls below, and freed afterwards.
  for (int i = 0; i < config->threads; i++) {
    svs.insert(svs.begin(), _p_svs[i].begin(), _p_svs[i].end());
    alignments.insert(alignments.begin(), _p_alignments[i].begin(),
                      _p_alignments[i].end());
  }
  sort(svs.begin(), svs.end());
  clean_dups();
  spdlog::info("{} SVs before chain filtering.", svs.size());
  // secondary sort by |svlen| so that same-position calls are adjacent by length;
  // this makes the chain filter deterministic and effective for same-position duplicates
  stable_sort(svs.begin(), svs.end(), [](const SV &a, const SV &b) {
    if (a.chrom != b.chrom) return false;
    if (a.s != b.s) return false;
    return abs(a.l) < abs(b.l);
  });
  filter_sv_chains();
  spdlog::info("Writing {} SVs.", svs.size());
  sort(svs.begin(), svs.end());
  write_vcf();

  if (config->poa.compare("") != 0) {
    spdlog::info("Writing POA alignments to {}..", config->poa);
    write_sam();
  }

  if (config->clipped) {
    spdlog::warn(
        "Calling imprecise SVs from clipped alignments is experimental");
    // Per-chromosome interval trees: ensures clips on chrom A cannot be
    // filtered/paired against SVs on chrom B that share a coordinate.
    unordered_map<string, interval_tree_t<int>> vartrees;
    for (const auto &sv : svs)
      vartrees[sv.chrom].insert({sv.s - 1000, sv.e + 1000});
    // Also exclude clips near germline-filtered events: these can produce
    // imprecise FP calls (e.g. reads soft-clipped at a germline insertion
    // boundary get paired as a spurious somatic INS or DEL).
    for (int i = 0; i < config->threads; i++)
      for (const auto &r : _p_germline_regions[i])
        vartrees[std::get<0>(r)].insert(
            {std::get<1>(r) - 1000, std::get<2>(r) + 1000});
    vector<SV> clipped_svs;
    Clipper clipper(C.clips);
    clipper.call(config->threads, vartrees);
    int s = 0;
    for (int i = 0; i < config->threads; i++) {
      s += clipper._p_svs[i].size();
      clipped_svs.insert(clipped_svs.begin(), clipper._p_svs[i].begin(),
                         clipper._p_svs[i].end());
    }
    spdlog::info("Predicted {} SVs from clipped alignments", s);
    // Deduplicate clipped SVs by reciprocal overlap (RO >= 0.9).
    // Left-clip and right-clip SA paths independently call the same large
    // deletion from opposite breakpoints; RO with max() denominator merges
    // them while preserving truly distinct events of different sizes.
    sort(clipped_svs.begin(), clipped_svs.end());
    vector<bool> suppressed(clipped_svs.size(), false);
    for (size_t i = 0; i < clipped_svs.size(); i++) {
      if (suppressed[i]) continue;
      for (size_t j = i + 1; j < clipped_svs.size(); j++) {
        if (suppressed[j]) continue;
        const SV &a = clipped_svs[i];
        const SV &b = clipped_svs[j];
        if (a.chrom != b.chrom || a.type != b.type) continue;
        int la = abs(a.l), lb = abs(b.l);
        int isect = max(0, min(a.s + la, b.s + lb) - max(a.s, b.s));
        double ro = (double)isect / max(la, lb);
        if (ro >= 0.9) {
          if (b.w > a.w) {
            suppressed[i] = true;
            break;
          } else {
            suppressed[j] = true;
          }
        }
      }
    }
    for (size_t i = 0; i < clipped_svs.size(); i++) {
      const SV &sv = clipped_svs[i];
      if (suppressed[i]) continue;
      // BND breakends have no meaningful length (l == 0); exempt them from the
      // min_sv_length filter, otherwise every translocation is silently dropped.
      if (sv.type == "BND" || abs(sv.l) >= (int)config->min_sv_length) {
        if (config->bed_filter.overlaps(sv.chrom, sv.s, sv.e)) continue;
        // Germline filter on every clipped event. Clipped calls come from a
        // read split (clip+SA), so the germline signal is a normal contig split
        // the same way — checked via the contigs' SA, not CIGAR I/D ops.
        // Serial loop → t=0.
        if (_p_normal_bam && is_germline_breakend(sv, 0))
          continue;
        cout << sv << endl;
      }
    }
  }

  if (_p_normal_bam) {
    for (int i = 0; i < config->threads; i++) {
      sam_hdr_destroy(_p_normal_hdr[i]);
      hts_idx_destroy(_p_normal_idx[i]);
      sam_close(_p_normal_bam[i]);
    }
    delete[] _p_normal_bam; _p_normal_bam = nullptr;
    delete[] _p_normal_idx; _p_normal_idx = nullptr;
    delete[] _p_normal_hdr; _p_normal_hdr = nullptr;
  }

  destroy_chromosomes();
}

void Caller::write_vcf() {
  print_vcf_header();
  for (const SV &sv : svs) {
    if (config->bed_filter.overlaps(sv.chrom, sv.s, sv.e)) {
      continue;
    }
    cout << sv << endl;
  }
}

void Caller::write_sam() {
  ofstream osam;
  osam.open(config->poa);
  osam << "@HD\tVN:1.4" << endl;
  for (size_t i = 0; i < chromosomes.size(); ++i)
    osam << "@SQ\tSN:" << chromosomes[i] << "\t"
         << "LN:" << strlen(chromosome_seqs[chromosomes[i]]) << endl;
  for (size_t j = 0; j < alignments.size(); j++)
    osam << alignments[j] << endl;
  osam.close();
}

// /* Split cluster in subclusters based on length */
vector<Cluster> Caller::split_cluster_by_len(const Cluster &cluster) {
  vector<Cluster> subclusters;
  for (uint c = 0; c < cluster.size(); ++c) {
    const SubRead &sr = cluster.get_subread(c);
    float sl = sr.size();
    size_t i;
    float ratio = -1.0f;
    for (i = 0; i < subclusters.size(); i++) {
      float cl = subclusters[i].get_len();
      float current_ratio = min(cl, sl) / max(cl, sl);
      ratio = current_ratio;
      if (current_ratio >= config->min_ratio)
        break;
    }
    if (i == subclusters.size()) {
      spdlog::debug(
          "[CALLER_FILTER][SPLIT_BY_LEN][NEW_SUBCLUSTER] chrom={} interval={} read={} len={} subclusters={} min_ratio={}",
          cluster.chrom, cluster.s, cluster.e, sr.name, sl, subclusters.size(),
          config->min_ratio);
      subclusters.push_back(Cluster(cluster.chrom, cluster.s, cluster.e,
                                    cluster.cov, cluster.cov0, cluster.cov1,
                                    cluster.cov2));
    } else {
      spdlog::debug(
          "[CALLER_FILTER][SPLIT_BY_LEN][ASSIGNED] chrom={} interval={} read={} len={} subcluster_index={} ratio={:.2f}",
          cluster.chrom, cluster.s, cluster.e, sr.name, sl, i, ratio);
    }
    subclusters[i].add_subread(sr);
  }
  return subclusters;
}

// Split cluster in subclusters
vector<Cluster> Caller::split_cluster(const Cluster &cluster) {
  // Step 1: split cluster by haplotype tag
  Cluster cluster_0 = cluster;
  cluster_0.clear();
  Cluster cluster_1 = cluster;
  cluster_1.clear();
  Cluster cluster_2 = cluster;
  cluster_2.clear();
  for (const SubRead &sr : cluster.subreads) {
    if (config->useht) {
      if (sr.htag == 1)
        cluster_1.add_subread(sr);
      else if (sr.htag == 2)
        cluster_2.add_subread(sr);
      else
        // 0 or no tag
        cluster_0.add_subread(sr);
    } else {
      cluster_0.add_subread(sr);
    }
  }
  cluster_0.cov1 = -1;
  cluster_0.cov2 = -1;
  cluster_1.cov0 = -1;
  cluster_1.cov2 = -1;
  cluster_2.cov0 = -1;
  cluster_2.cov1 = -1;

  spdlog::debug("[CALLER_SPLIT][HAP_DISTRIBUTION] chrom={} interval={} hap0={} hap1={} hap2={}",
                cluster.chrom, cluster.s, cluster.e, cluster_0.size(),
                cluster_1.size(), cluster_2.size());

  vector<Cluster> out_subclusters;
  if (cluster_1.size() == 0 && cluster_2.size() == 0) {
    // no alignment is tagged, use length
    vector<Cluster> subclusters = split_cluster_by_len(cluster_0);
    int i_max1 = -1, i_max2 = -1;
    uint v_max1 = 0, v_max2 = 0;
    for (uint i = 0; i < subclusters.size(); ++i) {
      if (subclusters[i].size() > v_max1) {
        v_max2 = v_max1;
        i_max2 = i_max1;
        v_max1 = subclusters[i].size();
        i_max1 = i;
      } else if (subclusters[i].size() > v_max2) {
        v_max2 = subclusters[i].size();
        i_max2 = i;
      }
    }
    string subcluster_sizes;
    for (const auto &sub : subclusters) {
      if (!subcluster_sizes.empty())
        subcluster_sizes += ",";
      subcluster_sizes += to_string(sub.size());
    }
    string kept1_desc = (i_max1 != -1)
                             ? "idx=" + to_string(i_max1) + " size=" +
                                   to_string(subclusters[i_max1].size())
                             : "none";
    string kept2_desc = (i_max2 != -1)
                             ? "idx=" + to_string(i_max2) + " size=" +
                                   to_string(subclusters[i_max2].size())
                             : "none";
    spdlog::debug(
        "[CALLER_SPLIT][LEN_ONLY] chrom={} interval={} subclusters={} sizes={} kept1={} kept2={}",
        cluster.chrom, cluster.s, cluster.e, subclusters.size(),
        subcluster_sizes, kept1_desc, kept2_desc);
    if (i_max1 != -1)
      out_subclusters.push_back(subclusters[i_max1]);
    if (i_max2 != -1)
      out_subclusters.push_back(subclusters[i_max2]);
  } else {
    int both = (cluster_1.size() > 0 ? 1 : 0) + (cluster_2.size() > 0 ? 2 : 0);
    vector<Cluster> subclusters_1 = split_cluster_by_len(cluster_1);
    vector<Cluster> subclusters_2 = split_cluster_by_len(cluster_2);
    Cluster new_cluster(cluster.chrom, cluster.s, cluster.e, cluster.cov,
                        cluster.cov0, -1, -1);

    for (uint c = 0; c < cluster_0.size(); ++c) {
      const SubRead &sr = cluster_0.get_subread(c);
      float sl = sr.size();

      int best_1 = -1;
      float best_ratio_1 = -1.0f;
      for (uint i = 0; i < subclusters_1.size(); i++) {
        float cl = subclusters_1[i].get_len();
        float r = min(cl, sl) / max(cl, sl);
        if (r >= config->min_ratio && r > best_ratio_1) {
          best_1 = i;
          best_ratio_1 = r;
        }
      }

      int best_2 = -1;
      float best_ratio_2 = -1.0f;
      for (uint i = 0; i < subclusters_2.size(); i++) {
        float cl = subclusters_2[i].get_len();
        float r = min(cl, sl) / max(cl, sl);
        if (r >= config->min_ratio && r > best_ratio_2) {
          best_2 = i;
          best_ratio_2 = r;
        }
      }

      if (both == 1) {
        assert(best_2 == -1);
        if (best_1 == -1) {
          new_cluster.add_subread(sr);
          spdlog::debug(
              "[CALLER_SPLIT][UNASSIGNED_HP1_ONLY] chrom={} interval={} read={} len={} new_cluster_reads={} ratio1={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl,
              new_cluster.size(), best_ratio_1);
        } else {
          subclusters_1[best_1].add_subread(sr);
          ++subclusters_1[best_1].cov1;
          --new_cluster.cov0;
          spdlog::debug(
              "[CALLER_SPLIT][ASSIGNED_HP1_ONLY] chrom={} interval={} read={} len={} hap1_subcluster={} ratio1={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl, best_1,
              best_ratio_1);
        }
      } else if (both == 2) {
        assert(best_1 == -1);
        if (best_2 == -1) {
          new_cluster.add_subread(sr);
          spdlog::debug(
              "[CALLER_SPLIT][UNASSIGNED_HP2_ONLY] chrom={} interval={} read={} len={} new_cluster_reads={} ratio2={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl,
              new_cluster.size(), best_ratio_2);
        } else {
          subclusters_2[best_2].add_subread(sr);
          ++subclusters_2[best_2].cov2;
          --new_cluster.cov0;
          spdlog::debug(
              "[CALLER_SPLIT][ASSIGNED_HP2_ONLY] chrom={} interval={} read={} len={} hap2_subcluster={} ratio2={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl, best_2,
              best_ratio_2);
        }
      } else {
        if (best_1 != -1 && best_ratio_1 > best_ratio_2) {
          subclusters_1[best_1].add_subread(sr);
          ++subclusters_1[best_1].cov1;
          --new_cluster.cov0;
          spdlog::debug(
              "[CALLER_SPLIT][ASSIGNED_BOTH] chrom={} interval={} read={} len={} hap1_subcluster={} ratio1={:.2f} ratio2={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl, best_1,
              best_ratio_1, best_ratio_2);
        } else if (best_2 != -1 && best_ratio_2 > best_ratio_1) {
          subclusters_2[best_2].add_subread(sr);
          ++subclusters_2[best_2].cov2;
          --new_cluster.cov0;
          spdlog::debug(
              "[CALLER_SPLIT][ASSIGNED_BOTH] chrom={} interval={} read={} len={} hap2_subcluster={} ratio1={:.2f} ratio2={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl, best_2,
              best_ratio_1, best_ratio_2);
        } else {
          new_cluster.add_subread(sr);
          spdlog::debug(
              "[CALLER_SPLIT][UNASSIGNED_BOTH] chrom={} interval={} read={} len={} new_cluster_reads={} ratio1={:.2f} ratio2={:.2f}",
              cluster.chrom, cluster.s, cluster.e, sr.name, sl,
              new_cluster.size(), best_ratio_1, best_ratio_2);
        }
      }
    }

    spdlog::debug(
        "[CALLER_SPLIT][UNASSIGNED_SUMMARY] chrom={} interval={} both={} shared_reads={} cov0={}",
        cluster.chrom, cluster.s, cluster.e, both, new_cluster.size(),
        new_cluster.cov0);

    uint v_max = 0;
    int i_max = -1;
    for (uint i = 0; i < subclusters_1.size(); ++i) {
      if (subclusters_1[i].size() > v_max) {
        v_max = subclusters_1[i].size();
        i_max = i;
      }
    }
    string hap1_sizes;
    for (const auto &sub : subclusters_1) {
      if (!hap1_sizes.empty())
        hap1_sizes += ",";
      hap1_sizes += to_string(sub.size());
    }
    string hap1_kept_desc = (i_max != -1)
                                  ? "idx=" + to_string(i_max) + " size=" +
                                        to_string(subclusters_1[i_max].size())
                                  : "none";
    spdlog::debug(
        "[CALLER_SPLIT][HAP1_FINAL] chrom={} interval={} hap1_subclusters={} sizes={} kept={}",
        cluster.chrom, cluster.s, cluster.e, subclusters_1.size(), hap1_sizes,
        hap1_kept_desc);
    if (i_max != -1)
      out_subclusters.push_back(subclusters_1[i_max]);
    v_max = 0, i_max = -1;
    for (uint i = 0; i < subclusters_2.size(); ++i) {
      if (subclusters_2[i].size() > v_max) {
        v_max = subclusters_2[i].size();
        i_max = i;
      }
    }
    string hap2_sizes;
    for (const auto &sub : subclusters_2) {
      if (!hap2_sizes.empty())
        hap2_sizes += ",";
      hap2_sizes += to_string(sub.size());
    }
    string hap2_kept_desc = (i_max != -1)
                                  ? "idx=" + to_string(i_max) + " size=" +
                                        to_string(subclusters_2[i_max].size())
                                  : "none";
    spdlog::debug(
        "[CALLER_SPLIT][HAP2_FINAL] chrom={} interval={} hap2_subclusters={} sizes={} kept={}",
        cluster.chrom, cluster.s, cluster.e, subclusters_2.size(), hap2_sizes,
        hap2_kept_desc);
    if (i_max != -1)
      out_subclusters.push_back(subclusters_2[i_max]);

    if (both != 3) {
      vector<Cluster> new_subclusters = split_cluster_by_len(new_cluster);
      string new_sizes;
      for (const auto &sub : new_subclusters) {
        if (!new_sizes.empty())
          new_sizes += ",";
        new_sizes += to_string(sub.size());
      }
      spdlog::debug(
          "[CALLER_SPLIT][UNASSIGNED_FALLBACK] chrom={} interval={} both={} new_cluster_reads={} new_subclusters={} sizes={}",
          cluster.chrom, cluster.s, cluster.e, both, new_cluster.size(),
          new_subclusters.size(), new_sizes);
      v_max = 0, i_max = -1;
      for (uint i = 0; i < new_subclusters.size(); ++i) {
        if (new_subclusters[i].size() > v_max) {
          v_max = new_subclusters[i].size();
          i_max = i;
        }
      }
      if (i_max != -1) {
        if (both == 1)
          new_subclusters[i_max].cov1 = -1;
        else
          new_subclusters[i_max].cov2 = -1;
        spdlog::debug(
            "[CALLER_SPLIT][UNASSIGNED_FALLBACK_KEEP] chrom={} interval={} selected_idx={} size={}",
            cluster.chrom, cluster.s, cluster.e, i_max,
            new_subclusters[i_max].size());
        out_subclusters.push_back(new_subclusters[i_max]);
      } else {
        spdlog::debug(
            "[CALLER_SPLIT][UNASSIGNED_FALLBACK_DROP] chrom={} interval={} no subcluster selected from fallback",
            cluster.chrom, cluster.s, cluster.e);
      }
    }
  }

  string result_sizes;
  for (const auto &sub : out_subclusters) {
    if (!result_sizes.empty())
      result_sizes += ",";
    result_sizes += to_string(sub.size());
  }
  spdlog::debug(
      "[CALLER_SPLIT][RESULT] chrom={} interval={} final_subclusters={} sizes={}",
      cluster.chrom, cluster.s, cluster.e, out_subclusters.size(), result_sizes);

  assert(out_subclusters.size() > 0 && out_subclusters.size() <= 2);
  return out_subclusters;
}

string Caller::run_poa(const vector<string> &seqs) {
  uint n_seqs = seqs.size();
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->align_mode = 0; // global
  abpt->disable_seeding = 1;
  abpt->progressive_poa = 0;
  abpt->amb_strand = 0;
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  // abpt->max_n_cons = 2;
  // abpt->min_freq = 0.25;
  abpoa_post_set_para(abpt);

  // abpt->match = 2;      // match score
  // abpt->mismatch = 4;   // mismatch penalty
  // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  // abpt->gap_open1 = 4;  // gap open penalty #1
  // abpt->gap_ext1 = 2;   // gap extension penalty #1
  // abpt->gap_open2 = 24; // gap open penalty #2
  // abpt->gap_ext2 = 1;   // gap extension penalty #2
  // gap_penalty = min{gap_open1 + gap_len*gap_ext1, gap_open2+gap_len*gap_ext2}

  int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
  uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
  for (uint i = 0; i < n_seqs; ++i) {
    seq_lens[i] = seqs[i].size();
    bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
    for (int j = 0; j < seq_lens[i]; ++j)
      bseqs[i][j] = _char26_table[(int)seqs[i][j]];
  }

  abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);
  abpoa_cons_t *abc = ab->abc;
  string cons = ""; // XXX: we may avoid converting to ACGT here since we need
                    // to reconvert back for ksw2
  if (abc->n_cons > 0)
    for (int j = 0; j < abc->cons_len[0]; ++j)
      cons += "ACGTN"[abc->cons_base[0][j]];

  for (uint i = 0; i < n_seqs; ++i)
    free(bseqs[i]);

  free(bseqs);
  free(seq_lens);
  abpoa_free(ab);
  abpoa_free_para(abpt);

  return cons;
}

// Returns true if sv is found in at least one normal contig that fully covers
// [cl_s, cl_e]. The contig's existing BAM CIGAR is parsed directly — no
// re-alignment needed since contigs are already mapped to the reference.
// All fully-covering contigs are checked to account for multiple alleles.
bool Caller::is_germline(const SV &sv, const string &chrom, int cl_s, int cl_e,
                         int t) {
  string region = chrom + ":" + to_string(cl_s) + "-" + to_string(cl_e);
  hts_itr_t *itr =
      sam_itr_querys(_p_normal_idx[t], _p_normal_hdr[t], region.c_str());
  if (!itr)
    return false;

  int sv_len = abs(sv.l);
  int sv_end = sv.s + sv_len;
  static const float min_ro = 0.95f;
  static const int min_len_diff = 2;
  bam1_t *aln = bam_init1();
  bool germline = false;

  while (!germline && sam_itr_next(_p_normal_bam[t], itr, aln) > 0) {
    if (aln->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
      continue;
    // Contig must fully cover the SV interval [sv.s, sv.s + sv_len].
    // Using the cluster boundaries [cl_s, cl_e] is too strict: the cluster
    // right boundary can be pulled far by unrelated SFS, excluding contigs
    // that perfectly span the SV but not the whole cluster.
    if (aln->core.pos > sv.s || (int)(bam_endpos(aln) - 1) < sv_end)
      continue;

    // Walk the contig's own CIGAR and check reciprocal interval overlap
    // between the called SV interval and each contig I/D operation.
    auto cigar_ops = decode_cigar(aln);
    uint rpos = aln->core.pos;
    for (const auto &op : cigar_ops) {
      uint l = op.first;
      int bam_op = op.second;
      if (bam_op == BAM_CMATCH || bam_op == BAM_CEQUAL ||
          bam_op == BAM_CDIFF) {
        rpos += l;
      } else if (bam_op == BAM_CINS || bam_op == BAM_CDEL) {
        bool type_match = (bam_op == BAM_CINS) ? (sv.type == "INS")
                                               : (sv.type == "DEL");
        if (type_match && (int)l >= (int)config->min_sv_length) {
          // SV interval: [sv.s, sv.s + sv_len]
          // contig op interval: [rpos, rpos + l]
          int overlap = max(0, min((int)sv.s + sv_len, (int)rpos + (int)l) -
                                  max((int)sv.s, (int)rpos));
          int max_len = max(sv_len, (int)l);
          float ro = (max_len > 0) ? (float)overlap / max_len : 0.0f;
          spdlog::info(
              "[GERMLINE_CHECK] sv={}:{} type={} sv_len={} | "
              "contig={} op={} op_rpos={} op_len={} | overlap={} ro={:.4f} threshold={:.4f} -> {}",
              chrom, sv.s, sv.type, sv_len,
              bam_get_qname(aln), (bam_op == BAM_CINS ? 'I' : 'D'), rpos, l,
              overlap, ro, min_ro, (ro >= min_ro ? "GERMLINE" : "no_match"));
          if (ro >= min_ro || (abs((int)l - sv_len) <= min_len_diff && ro > 0.9f)) {
            germline = true;
            break;
          }
        }
        if (bam_op == BAM_CDEL)
          rpos += l;
        // early exit once we are clearly past the SV region
        if ((int)rpos > sv.s + sv_len + sv_len)
          break;
      }
    }
  }

  bam_destroy1(aln);
  hts_itr_destroy(itr);
  return germline;
}

// Germline check for clipped-SFS calls (any type). SVDSS called these from a
// read clip whose SA defines the event; the germline equivalent is a normal
// CONTIG split the same way: clipped near the breakpoint with an SA reaching the
// same partner. We confirm germline if some normal contig near sv.s has an SA
// congruent with the SA SVDSS chose. The partner and the required strand
// relationship are reconstructed from the SV geometry per type:
//   - BND: mate chrom:pos from the breakend ALT ("N]chr5:1234]"); any strand.
//   - INV: same chrom, sv.s + |l|; opposite strand (inverted junction).
//   - DEL/DUP: same chrom, sv.s + |l|; same strand.
//   - INS: same chrom, ~sv.s; same strand.
bool Caller::is_germline_breakend(const SV &sv, int t) {
  if (!_p_normal_bam || !_p_normal_idx[t] || !_p_normal_hdr[t])
    return false;

  string mate_chrom;
  long mate_pos = -1;
  // required strand relationship of the contig split: 0 = any, 1 = same, 2 = flip
  int strand_req = 1;
  if (sv.type == "BND") {
    // parse "chrom:pos" inside the bracket notation of altall
    size_t lb = sv.altall.find_first_of("[]");
    size_t rb = sv.altall.find_first_of("[]", lb + 1);
    if (lb == string::npos || rb == string::npos)
      return false;
    string inside = sv.altall.substr(lb + 1, rb - lb - 1); // chrom:pos
    size_t colon = inside.rfind(':');
    if (colon == string::npos)
      return false;
    mate_chrom = inside.substr(0, colon);
    try {
      mate_pos = stol(inside.substr(colon + 1));
    } catch (...) {
      return false;
    }
    strand_req = 0; // translocation orientation not constrained
  } else if (sv.type == "INV") {
    mate_chrom = sv.chrom;
    mate_pos = (long)sv.s + abs(sv.l);
    strand_req = 2;
  } else if (sv.type == "DEL" || sv.type == "DUP") {
    mate_chrom = sv.chrom;
    mate_pos = (long)sv.s + abs(sv.l);
    strand_req = 1;
  } else { // INS
    mate_chrom = sv.chrom;
    mate_pos = sv.s;
    strand_req = 1;
  }
  if (mate_pos < 0)
    return false;

  const int W = 500; // position tolerance (contig vs read, HPC, breakpoint)
  int qs = sv.s > W ? sv.s - W : 0;
  string region =
      sv.chrom + ":" + to_string(qs) + "-" + to_string(sv.s + W);
  hts_itr_t *itr =
      sam_itr_querys(_p_normal_idx[t], _p_normal_hdr[t], region.c_str());
  if (!itr)
    return false;

  bam1_t *aln = bam_init1();
  bool germline = false;
  while (!germline && sam_itr_next(_p_normal_bam[t], itr, aln) > 0) {
    if (aln->core.flag & BAM_FSECONDARY)
      continue;
    uint8_t *sa_tag = bam_aux_get(aln, "SA");
    if (sa_tag == NULL)
      continue;
    bool contig_reverse = (aln->core.flag & BAM_FREVERSE) != 0;
    // Each SA entry: "chrom,pos,strand,CIGAR,mapQ,NM;"
    string sa(bam_aux2Z(sa_tag));
    size_t s = 0;
    while (!germline && s < sa.size()) {
      size_t e = sa.find(';', s);
      if (e == string::npos)
        break;
      string entry = sa.substr(s, e - s);
      s = e + 1;
      // split the 6 comma fields
      size_t c1 = entry.find(','), c2 = entry.find(',', c1 + 1),
             c3 = entry.find(',', c2 + 1), c4 = entry.find(',', c3 + 1);
      if (c1 == string::npos || c2 == string::npos || c3 == string::npos ||
          c4 == string::npos)
        continue;
      string sa_chrom = entry.substr(0, c1);
      long sa_pos;
      try {
        sa_pos = stol(entry.substr(c1 + 1, c2 - c1 - 1));
      } catch (...) {
        continue;
      }
      bool sa_reverse = (entry[c2 + 1] == '-');
      // reference span of the SA from its CIGAR (M/=/X/D/N consume reference)
      string sa_cigar = entry.substr(c3 + 1, c4 - c3 - 1);
      long sa_ref_len = 0, num = 0;
      for (char ch : sa_cigar) {
        if (isdigit((unsigned char)ch))
          num = num * 10 + (ch - '0');
        else {
          if (ch == 'M' || ch == '=' || ch == 'X' || ch == 'D' || ch == 'N')
            sa_ref_len += num;
          num = 0;
        }
      }
      if (sa_chrom != mate_chrom)
        continue;
      bool flip = (sa_reverse != contig_reverse);
      if (strand_req == 2 && !flip) // need an inverted junction
        continue;
      if (strand_req == 1 && flip) // need a same-strand junction
        continue;
      // does the SA alignment reach the partner breakpoint?
      if (sa_pos - W <= mate_pos && mate_pos <= sa_pos + sa_ref_len + W)
        germline = true;
    }
  }
  bam_destroy1(aln);
  hts_itr_destroy(itr);
  return germline;
}

// Call SVs by POA+realignment
void Caller::pcall(const vector<Cluster> &clusters) {
  // Aggregate statistics (info-level summary at the end of the step).
  size_t n_bed_excluded = 0;   // clusters dropped by the BED filter
  size_t n_below_weight = 0;   // clusters with size < min_cluster_weight
  size_t n_passed_weight = 0;  // clusters surviving both filters (eligible)
  size_t n_with_call = 0;      // clusters that produced >=1 reported SV
  // Sub-cluster level breakdown of what happens to eligible clusters.
  size_t n_subclusters = 0;        // sub-clusters produced by the split
  size_t n_sub_below_weight = 0;   // sub-clusters below threshold after split
  size_t n_sub_no_sv = 0;          // sub-clusters whose POA yielded no SV >=min_sv_length
  size_t n_cand_svs = 0;           // candidate SVs from POA (before SV-level filters)
  size_t n_filt_lowoverlap = 0;    // candidate SVs dropped: low original-SFS overlap
  size_t n_filt_germline = 0;      // candidate SVs dropped: germline match
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < clusters.size(); i++) {
    int t = omp_get_thread_num();
    const Cluster &cluster = clusters[i];

    if (config->bed_filter.overlaps(cluster.chrom, cluster.s, cluster.e)) {
      spdlog::debug(
          "[CALLER_FILTER][BED_EXCLUSION] chrom={} interval={}-{} reads={}",
          cluster.chrom, cluster.s, cluster.e, cluster.size());
#pragma omp atomic
      ++n_bed_excluded;
      continue;
    }

    if (cluster.size() < config->min_cluster_weight) {
      spdlog::debug(
          "[CALLER_FILTER][MIN_CLUSTER_WEIGHT] chrom={} interval={} reads={} threshold={}",
          cluster.chrom, cluster.s, cluster.e, cluster.size(),
          config->min_cluster_weight);
#pragma omp atomic
      ++n_below_weight;
      continue;
    }
#pragma omp atomic
    ++n_passed_weight;
    size_t svs_before = _p_svs[t].size(); // to detect if this cluster calls
    string chrom = cluster.chrom;

    const vector<Cluster> &subclusters = split_cluster(cluster);
    string split_sizes;
    for (const auto &sub : subclusters) {
      if (!split_sizes.empty())
        split_sizes += ",";
      split_sizes += to_string(sub.size());
    }
    spdlog::debug(
        "[CALLER_SPLIT][RESULT] chrom={} interval={}-{} original_reads={} subclusters={} sizes={}",
        cluster.chrom, cluster.s, cluster.e, cluster.size(), subclusters.size(),
        split_sizes);
#pragma omp atomic
    n_subclusters += subclusters.size();

    // Calling from one or two clusters
    for (const Cluster &cl : subclusters) {
      if (cl.size() < config->min_cluster_weight) {
        spdlog::debug(
            "[CALLER_FILTER][MIN_CLUSTER_WEIGHT_POST_SPLIT] chrom={} interval={} subcluster_reads={} threshold={}",
            cluster.chrom, cluster.s, cluster.e, cl.size(),
            config->min_cluster_weight);
#pragma omp atomic
        ++n_sub_below_weight;
      }

      vector<SV> _svs;

      string ref = string(chromosome_seqs[chrom] + cl.s, cl.e - cl.s + 1);
      string consensus = run_poa(cl.get_seqs());

      // ksw2 stuff - TODO: move to a separate function
      int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
      int8_t a = (int8_t)sc_mch,
             b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
      int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                        b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
      uint tl = ref.size(), ql = consensus.size();
      uint8_t *ts = (uint8_t *)malloc(tl);
      uint8_t *qs = (uint8_t *)malloc(ql);
      for (i = 0; i < tl; ++i)
        ts[i] = _char26_table[(uint8_t)ref[i]]; // encode to 0/1/2/3
      for (i = 0; i < ql; ++i)
        qs[i] = _char26_table[(uint8_t)consensus[i]];

      ksw_extz_t ez;
      memset(&ez, 0, sizeof(ksw_extz_t));
      ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1,
                    -1, 0, &ez);

      int score = ez.score;
      string cigar_str = "";
      for (int i = 0; i < ez.n_cigar; ++i) {
        cigar_str += to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
      }
      _p_alignments[t].push_back(
          Consensus(consensus, cigar_str, chrom, cl.s, cl.e));

      // -- Extracting SVs
      uint rpos = cl.s; // position on reference
      uint cpos = 0;    // position on consensus
      CIGAR cigar;
      cigar.parse_cigar(cigar_str.c_str());
      int nv = 0;
      for (const auto &cigar_pair : cigar.ops) {
        uint l = cigar_pair.first;
        char op = cigar_pair.second;
        if (op == '=' || op == 'M') {
          rpos += l;
          cpos += l;
        } else if (op == 'I') {
          if (l >= config->min_sv_length) {
            SV sv = SV("INS", cl.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, 1),
                       string(chromosome_seqs[chrom] + rpos - 1, 1) +
                           consensus.substr(cpos, l),
                       cl.size(), cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(cl.get_names());
            _svs.push_back(sv);
            nv++;
          }
          cpos += l;
        } else if (op == 'D') {
          if (l >= config->min_sv_length) {
            SV sv = SV("DEL", cl.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, l),
                       string(chromosome_seqs[chrom] + rpos - 1, 1), cl.size(),
                       cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(cl.get_names());
            _svs.push_back(sv);
            nv++;
          }
          rpos += l;
        }
      }
      if (_svs.empty()) {
#pragma omp atomic
        ++n_sub_no_sv;
      } else {
#pragma omp atomic
        n_cand_svs += _svs.size();
      }
      for (size_t v = 0; v < _svs.size(); v++) {
        _svs[v].ngaps = nv;
        _svs[v].set_gt("./.", 100);
        _svs[v].set_cov(cl.cov, cl.cov0, cl.cov1, cl.cov2);
        _svs[v].set_rvec(cluster.reads);
      }
      for (const SV &sv : _svs) {
        if (config->require_sfs_overlap &&
            count_overlapping_original_sfs(cluster, cl, sv) <
                2) {
          spdlog::debug(
              "[CALLER_FILTER][LOW_ORIG_SFS_OVERLAP] chrom={} sv_start={} sv_end={} type={} cluster_interval={}:{}-{}",
              sv.chrom, sv.s, sv.e, sv.type, cluster.chrom, cluster.s,
              cluster.e);
#pragma omp atomic
          ++n_filt_lowoverlap;
          continue;
        }
        if (_p_normal_bam && is_germline(sv, chrom, cl.s, cl.e, t)) {
          spdlog::debug(
              "[CALLER_FILTER][GERMLINE] chrom={} sv_start={} sv_end={} type={} len={}",
              sv.chrom, sv.s, sv.e, sv.type, sv.l);
#pragma omp atomic
          ++n_filt_germline;
          // Record the germline region so that clip-based calls near this
          // event are also suppressed by filter_tooclose_clips().
          _p_germline_regions[t].push_back(
              std::make_tuple(sv.chrom, sv.s, sv.s + abs(sv.l)));
          continue;
        }
        _p_svs[t].push_back(sv);
      }
    }
    if (_p_svs[t].size() > svs_before) {
#pragma omp atomic
      ++n_with_call;
    }
  }
  spdlog::info(
      "[CALLER_STATS] clusters: total={} bed_excluded={} below_weight(<{})={} "
      "passed_weight={} produced_call={}",
      clusters.size(), n_bed_excluded, config->min_cluster_weight,
      n_below_weight, n_passed_weight, n_with_call);
  spdlog::info(
      "[CALLER_STATS] subclusters: total={} below_weight_post_split={} "
      "no_sv_from_poa={} | candidate_svs={} filtered_lowoverlap={} "
      "filtered_germline={} kept={}",
      n_subclusters, n_sub_below_weight, n_sub_no_sv, n_cand_svs,
      n_filt_lowoverlap, n_filt_germline,
      n_cand_svs - n_filt_lowoverlap - n_filt_germline);
}

int Caller::count_overlapping_original_sfs(const Cluster &cluster,
                                           const Cluster &subcluster,
                                           const SV &sv) const {
  const vector<string> sub_reads = subcluster.get_names();
  if (sub_reads.empty())
    return 0;
  unordered_set<string> subread_set(sub_reads.begin(), sub_reads.end());
  int count = 0;
  for (const SFS &sfs : cluster.SFSs) {
    if (subread_set.find(sfs.qname) == subread_set.end())
      continue;
    bool overlaps = false;
    if (!sfs.orig_intervals.empty()) {
      for (const auto &interval : sfs.orig_intervals) {
        if (interval.first <= sv.e && interval.second >= sv.s) {
          overlaps = true;
          break;
        }
      }
    } else {
      if (sfs.rs <= sv.e && sfs.re >= sv.s)
        overlaps = true;
    }
    if (overlaps)
      count++;
  }
  return count;
}

// Clean same SV reported twice
void Caller::clean_dups() {
  vector<SV> _svs;
  string last_chrom = "";
  int last_pos = -1;
  string last_refall = "";
  string last_altall = "";
  for (size_t i = 0; i < svs.size(); i++) {
    if (last_chrom != svs[i].chrom || last_pos != svs[i].s ||
        last_refall != svs[i].refall || last_altall != svs[i].altall)
      _svs.push_back(svs[i]);
    last_chrom = svs[i].chrom;
    last_pos = svs[i].s;
    last_refall = svs[i].refall;
    last_altall = svs[i].altall;
  }
  svs.clear();
  svs.insert(svs.begin(), _svs.begin(), _svs.end());
}

/* Merge close and similar SVs */
void Caller::filter_sv_chains() {
  if (svs.size() < 2)
    return;

  vector<SV> _svs;
  SV prev = svs[0];
  for (size_t i = 1; i < svs.size(); i++) {
    const SV &sv = svs[i];
    bool merged = false;
    if (sv.chrom == prev.chrom && sv.s - prev.e < 2 * sv.l &&
        prev.type == sv.type) {
      double w_r =
          min((double)sv.w, (double)prev.w) / max((double)sv.w, (double)prev.w);
      double l_r =
          min((double)sv.l, (double)prev.l) / max((double)sv.l, (double)prev.l);
      int d = sv.s - prev.s;
      // bypass weight ratio check when calls land at the same position:
      // they are duplicates from different clusters of the same event
      bool weight_ok = (d == 0) || (w_r >= 0.9);
      if (d < 100 && weight_ok &&
          l_r >= config->min_ratio) {
        double sim;
        if (sv.type == "DEL")
          sim = rapidfuzz::fuzz::ratio(sv.refall, prev.refall);
        else
          sim = rapidfuzz::fuzz::ratio(sv.altall, prev.altall);
        if (sim > 70) {
          spdlog::debug("[CALLER_CHAIN][MERGE] chrom={} prev_start={} prev_w={} new_start={} new_w={} type={} sim={:.1f} keep={}",
                        sv.chrom, prev.s, prev.w, sv.s, sv.w, sv.type, sim,
                        sv.w > prev.w ? "new" : "prev");
          // keep the winner as prev and continue comparing — this correctly
          // handles chains of 3+ duplicates without the reset-skip bug
          if (sv.w > prev.w)
            prev = sv;
          merged = true;
        }
      }
    }
    if (!merged) {
      _svs.push_back(prev);
      prev = sv;
    }
  }
  _svs.push_back(prev);
  svs.clear();
  svs.insert(svs.begin(), _svs.begin(), _svs.end());
}

void Caller::print_vcf_header() {
  cout << "##fileformat=VCFv4.2" << endl;
  cout << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
          "data_collections/HGSVC2/technical/reference/"
          "20200513_hg38_NoALT/"
          "hg38.no_alt.fa.gz"
       << endl;
  for (size_t i = 0; i < chromosomes.size(); ++i) {
    cout << "##contig=<ID=" << chromosomes[i]
         << ",length=" << strlen(chromosome_seqs[chromosomes[i]]) << ">"
         << endl;
  }
  cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
  cout << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant "
          "class\">"
       << endl;
  cout << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant "
          "type\">"
       << endl;
  cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="
          "\"Difference in "
          "length between REF and ALT alleles\">"
       << endl;
  cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End "
          "position of "
          "the variant described in this record\">"
       << endl;
  cout << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number "
          "of "
          "alignments supporting this record\">"
       << endl;
  cout << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus\">"
       << endl;
  cout << "##INFO=<ID=COV0,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (no HP)\">"
       << endl;
  cout << "##INFO=<ID=COV1,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (HP=1)\">"
       << endl;
  cout << "##INFO=<ID=COV2,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (HP=2)\">"
       << endl;
  cout << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Alignment "
          "score\">"
       << endl;
  cout << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of "
          "variations on same consensus\">"
       << endl;
  cout << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="
          "\"Imprecise "
          "structural variation\">"
       << endl;
  cout << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR of "
          "consensus\">"
       << endl;
  cout << "##INFO=<ID=READS,Number=.,Type=String,Description=\"Reads "
          "identifiers supporting the call\">"
       << endl;
    cout << "##INFO=<ID=SA_READS,Number=.,Type=String,Description=\"Reads "
      "with supplementary alignments supporting the call\">"
         << endl;
  cout << "##INFO=<ID=RVEC,Number=.,Type=String,Description=\"Reads vector "
          "used by genotyper\">"
       << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       << endl;
  cout << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype "
          "quality\">"
       << endl;
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT"
       << endl;
}
