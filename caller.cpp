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
  } else {
    // Every germline knob (min-ro, min-reads, max-reads, diff, realign) reads
    // the normal BAM, so without --normal-contigs-bam the whole filter is a
    // no-op. Say so: passing the flags and silently not filtering is worse than
    // not passing them at all.
    spdlog::warn("No --normal-contigs-bam: the germline filter is DISABLED and "
                 "every --germline-* flag is ignored. All calls are reported as "
                 "somatic.");
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
        if (excluded_by_bed_or_N(sv)) continue;
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
    if (excluded_by_bed_or_N(sv)) {
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

vector<PoaCons> Caller::run_poa(const vector<string> &seqs) {
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
  // Diploid mode: let abPOA emit up to two consensus sequences (one per allele)
  // when a heterozygous variant exceeds poa_min_freq. Splits blended VNTR/multi-
  // allele clusters that a single consensus averages into a phantom size.
  if (config->diploid_poa) {
    abpt->max_n_cons = 2;
    abpt->min_freq = config->poa_min_freq;
  }
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

  vector<PoaCons> out;
  for (int c = 0; c < abc->n_cons; ++c) {
    PoaCons pc;
    for (int j = 0; j < abc->cons_len[c]; ++j)
      pc.seq += "ACGTN"[abc->cons_base[c][j]];
    // Read membership: abPOA fills clu_read_ids only when it emits >1 consensus;
    // with a single consensus every read belongs to it.
    if (abc->n_cons > 1 && abc->clu_n_seq && abc->clu_read_ids) {
      for (int r = 0; r < abc->clu_n_seq[c]; ++r)
        pc.read_ids.push_back(abc->clu_read_ids[c][r]);
    } else {
      for (uint r = 0; r < n_seqs; ++r)
        pc.read_ids.push_back((int)r);
    }
    out.push_back(std::move(pc));
  }

  for (uint i = 0; i < n_seqs; ++i)
    free(bseqs[i]);
  free(bseqs);
  free(seq_lens);
  abpoa_free(ab);
  abpoa_free_para(abpt);

  // Safeguard: two consensuses that differ in length by less than min_sv_length
  // cannot yield distinct SVs (the size difference is below what we'd call). We
  // have no evidence to separate them, so merge into the higher-support one and
  // pool the reads (full weight). Prevents phantom pairs from spurious 1-2 bp
  // abPOA splits of unimodal clusters (e.g. the 120-vs-121 bp case).
  if (out.size() == 2 &&
      abs((int)out[0].seq.size() - (int)out[1].seq.size()) <
          (int)config->min_sv_length) {
    size_t maj = out[0].read_ids.size() >= out[1].read_ids.size() ? 0 : 1;
    size_t mnr = 1 - maj;
    out[maj].read_ids.insert(out[maj].read_ids.end(), out[mnr].read_ids.begin(),
                             out[mnr].read_ids.end());
    PoaCons keep = std::move(out[maj]);
    out.clear();
    out.push_back(std::move(keep));
  }
  return out;
}

// Align query to ref with the SAME ksw2 scoring/mode as the tumor consensus
// (global, asm20-like: match 1, mis -4, gap 6/26, ext 2/1). Returns the CIGAR as
// (length, op) pairs with op 0=M, 1=I, 2=D. Re-aligning normal reads through this
// makes their repeat gap-fragmentation match the tumor side.
vector<pair<int, int>> Caller::ksw_align(const string &ref, const string &query) {
  vector<pair<int, int>> ops;
  uint tl = ref.size(), ql = query.size();
  if (tl == 0 || ql == 0)
    return ops;
  int gapo = 6, gape = 2, gapo2 = 26, gape2 = 1;
  int8_t a = 1, b = -4;
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  uint8_t *ts = (uint8_t *)malloc(tl);
  uint8_t *qs = (uint8_t *)malloc(ql);
  for (uint i = 0; i < tl; ++i)
    ts[i] = _char26_table[(uint8_t)ref[i]];
  for (uint i = 0; i < ql; ++i)
    qs[i] = _char26_table[(uint8_t)query[i]];
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1, -1,
                0, &ez);
  for (int i = 0; i < ez.n_cigar; ++i)
    ops.push_back({(int)(ez.cigar[i] >> 4), (int)(ez.cigar[i] & 0xf)});
  free(ts);
  free(qs);
  free(ez.cigar);
  return ops;
}

// Re-align a single normal read's window around sv with ksw_align and report
// whether it carries an indel of sv's type with concordant size. This defeats the
// minimap2-vs-ksw2 gap-placement mismatch that lets fragmented repeat calls evade
// the germline filter.
bool Caller::normal_read_concordant(bam1_t *aln, const string &chrom,
                                    const SV &sv, int cl_s, int cl_e, int want,
                                    float min_ratio_len, int diff_max,
                                    int diff_min_len) {
  auto it = chromosome_seqs.find(chrom);
  if (it == chromosome_seqs.end() || it->second == nullptr)
    return false;
  const char *cseq = it->second;
  long clen = (long)strlen(cseq);
  int sv_len = abs(sv.l);
  // Align over the SAME window the tumor consensus used ([cl_s, cl_e]) so the
  // gap fragmentation in repeats is reproduced identically. Sizing the window to
  // the cluster (not a fixed margin) keeps it tight where clusters are small (so
  // most reads still span it) and wide only in the long repeats that fragment.
  long win_s = (long)cl_s;
  if (win_s < 0)
    win_s = 0;
  long win_e = (long)cl_e + 1;
  if (win_e > clen)
    win_e = clen;
  if (win_e <= win_s)
    return false;
  // read must fully span the window for a clean local re-alignment
  if ((long)aln->core.pos > win_s || (long)(bam_endpos(aln) - 1) < win_e - 1)
    return false;

  // extract the read bases aligned within [win_s, win_e)
  auto cigar_ops = decode_cigar(aln);
  uint8_t *qseq = bam_get_seq(aln);
  long rpos = aln->core.pos; // ref (0-based)
  long qpos = 0;             // query index
  string read_win;
  for (const auto &op : cigar_ops) {
    int l = op.first, bam_op = op.second;
    if (bam_op == BAM_CMATCH || bam_op == BAM_CEQUAL || bam_op == BAM_CDIFF) {
      for (int k = 0; k < l; ++k) {
        if (rpos >= win_s && rpos < win_e)
          read_win += seq_nt16_str[bam_seqi(qseq, qpos)];
        rpos++;
        qpos++;
      }
    } else if (bam_op == BAM_CINS) {
      if (rpos > win_s && rpos <= win_e)
        for (int k = 0; k < l; ++k)
          read_win += seq_nt16_str[bam_seqi(qseq, qpos + k)];
      qpos += l;
    } else if (bam_op == BAM_CDEL || bam_op == BAM_CREF_SKIP) {
      rpos += l;
    } else if (bam_op == BAM_CSOFT_CLIP) {
      qpos += l;
    } // BAM_CHARD_CLIP/PAD: consume nothing
    if (rpos >= win_e)
      break;
  }
  if (read_win.empty())
    return false;

  string ref_win(cseq + win_s, win_e - win_s);
  auto ops2 = ksw_align(ref_win, read_win);

  long rp = win_s; // absolute ref position while walking the re-aligned CIGAR
  for (const auto &o : ops2) {
    int l = o.first, op = o.second; // op: 0=M, 1=I(=BAM_CINS), 2=D(=BAM_CDEL)
    if (op == want) {
      const int nl = l;
      const bool ratio_ok =
          nl >= (int)config->min_sv_length &&
          (float)min(nl, sv_len) / (float)max(nl, sv_len) >= min_ratio_len;
      const bool diff_ok = diff_max > 0 && nl >= diff_min_len &&
                           abs(sv_len - nl) < diff_max;
      if (ratio_ok || diff_ok)
        return true;
    }
    if (op == 0 || op == BAM_CDEL) // M or D consume reference
      rp += l;
  }
  return false;
}

// Read-based germline filter. An SV (INS/DEL) is called germline when at least
// config.germline_min_reads normal reads carry a concordant indel near the
// breakpoint: same type and length ratio >= config.germline_min_ro. Only reads
// with MAPQ >= config.min_mapq (the SVDSS-call threshold) are considered, and at
// most config.germline_max_reads reads are examined per SV to bound cost in
// high-depth/repeat loci (early-exit once enough support is found). The normal
// BAM (config.normal_contigs_bam) may be a normal *reads* BAM; unlike the old
// contig-CIGAR check this sees per-haplotype/per-allele read evidence that the
// assembled contigs miss (e.g. VNTR alleles absent from the unitig set).
bool Caller::is_germline(const SV &sv, const string &chrom, int cl_s, int cl_e,
                         int t) {
  const int want = (sv.type == "INS") ? BAM_CINS
                   : (sv.type == "DEL") ? BAM_CDEL
                                        : -1;
  if (want == -1)
    return false; // read-based filter only handles INS/DEL

  int sv_len = abs(sv.l);
  int sv_end = sv.s + sv_len;
  const int win = 150; // window around the breakpoint
  string region = chrom + ":" + to_string(max(0, sv.s - win)) + "-" +
                  to_string(sv_end + win);
  hts_itr_t *itr =
      sam_itr_querys(_p_normal_idx[t], _p_normal_hdr[t], region.c_str());
  if (!itr)
    return false;

  const float min_ratio_len = config->germline_min_ro;
  // Difference-based germline test (union with the ratio test). diff_max of 0
  // means "use min_sv_length" (so it scales with the calling threshold).
  const int diff_max = config->germline_diff
                           ? (config->germline_diff_max > 0
                                  ? config->germline_diff_max
                                  : (int)config->min_sv_length)
                           : 0;
  const int diff_min_len = config->germline_diff_min_len;
  bam1_t *aln = bam_init1();
  int support = 0;  // normal reads carrying a concordant indel
  int examined = 0; // normal reads passing the flag/MAPQ filter
  bool germline = false;

  while (!germline && sam_itr_next(_p_normal_bam[t], itr, aln) > 0) {
    if (aln->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
      continue;
    if ((int)aln->core.qual < (int)config->min_mapq)
      continue; // only reads with at least the SVDSS-call min mapq
    if (++examined > config->germline_max_reads)
      break; // cap: bound cost in high-depth/repeat loci
    // read must span the breakpoint to carry the event
    if (aln->core.pos > sv.s || (int)(bam_endpos(aln) - 1) < sv.s)
      continue;

    // A read is concordant if its BAM (minimap2) CIGAR carries the indel OR, with
    // --germline-realign, re-aligning its window with the tumor's ksw2 does. The
    // realign is ADDITIVE (only a fallback when the BAM scan fails), so it can only
    // ADD detections (repeat calls the tumor fragments and the map-hifi CIGAR
    // misses) and never filters fewer reads than the baseline.
    bool read_concordant = false;
    auto cigar_ops = decode_cigar(aln);
    uint rpos = aln->core.pos;
    for (const auto &op : cigar_ops) {
      uint l = op.first;
      int bam_op = op.second;
      bool consumes_ref =
          (bam_op == BAM_CMATCH || bam_op == BAM_CEQUAL ||
           bam_op == BAM_CDIFF || bam_op == BAM_CDEL || bam_op == BAM_CREF_SKIP);
      if (bam_op == want && (int)rpos >= sv.s - win &&
          (int)rpos <= sv_end + win) {
        const int nl = (int)l;
        // Ratio test: normal indel of comparable relative length (>= min_sv_length).
        const bool ratio_ok =
            nl >= (int)config->min_sv_length &&
            (float)min(nl, sv_len) / (float)max(nl, sv_len) >= min_ratio_len;
        // Difference test: normal indel present but a different number of repeat
        // units; catches VNTR/low-complexity jitter the ratio test misses.
        const bool diff_ok = diff_max > 0 && nl >= diff_min_len &&
                             abs(sv_len - nl) < diff_max;
        if (ratio_ok || diff_ok) {
          read_concordant = true;
          break; // one concordant indel per read is enough
        }
      }
      if (consumes_ref)
        rpos += l;
      if ((int)rpos > sv_end + win)
        break; // past the region
    }
    if (!read_concordant && config->germline_realign &&
        (cl_e - cl_s) <= config->germline_realign_max_len)
      read_concordant = normal_read_concordant(
          aln, chrom, sv, cl_s, cl_e, want, min_ratio_len, diff_max, diff_min_len);
    if (read_concordant)
      ++support;
    if (support >= config->germline_min_reads)
      germline = true;
  }

  spdlog::debug("[GERMLINE_CHECK] sv={}:{} type={} sv_len={} | normal_examined={} "
                "concordant={} (min_reads={} min_ratio={:.2f} mapq>={}) -> {}",
                chrom, sv.s, sv.type, sv_len, examined, support,
                config->germline_min_reads, min_ratio_len, config->min_mapq,
                germline ? "GERMLINE" : "somatic");

  bam_destroy1(aln);
  hts_itr_destroy(itr);
  return germline;
}

// True if the reference has an N within +-window bp of a (1-based) position.
// chromosome_seqs holds 0-based, null-terminated sequences; a missing chrom
// (e.g. an ALT-contig mate absent from the reference) yields false.
bool Caller::ref_has_N_near(const string &chrom, long pos, int window) const {
  auto it = chromosome_seqs.find(chrom);
  if (it == chromosome_seqs.end() || it->second == nullptr)
    return false;
  const char *seq = it->second;
  long len = (long)strlen(seq);
  if (len == 0)
    return false;
  long c = pos - 1; // 1-based VCF pos -> 0-based index
  long s = c - window; if (s < 0) s = 0;
  long e = c + window; if (e > len - 1) e = len - 1;
  for (long i = s; i <= e; i++)
    if (seq[i] == 'N' || seq[i] == 'n')
      return true;
  return false;
}

// Parse the mate chrom:pos out of a BND ALT such as "N]chr5:1234]" or
// "[chr5:1234[T". Returns false for non-BND or unparseable ALTs.
bool Caller::bnd_mate(const SV &sv, string &mate_chrom, long &mate_pos) const {
  if (sv.type != "BND")
    return false;
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
  return true;
}

// Central exclusion test for a called SV. Drops it when either breakpoint falls
// in the exclusion bed OR lies within +-100 bp of a reference N gap. For BND
// both the primary and the mate breakpoint are tested (the bed check historically
// looked only at the primary end, letting subtelomeric cross-mapping BNDs whose
// mate lands in a gap slip through).
bool Caller::excluded_by_bed_or_N(const SV &sv) const {
  const int NWIN = 100;
  // primary breakpoint(s): sv.s and sv.e cover both ends of DEL/DUP/INV.
  if (config->bed_filter.overlaps(sv.chrom, sv.s, sv.e))
    return true;
  if (ref_has_N_near(sv.chrom, sv.s, NWIN) ||
      ref_has_N_near(sv.chrom, sv.e, NWIN))
    return true;
  if (sv.type == "BND") {
    string mc;
    long mp;
    if (bnd_mate(sv, mc, mp)) {
      if (config->bed_filter.overlaps(mc, mp, mp))
        return true;
      if (ref_has_N_near(mc, mp, NWIN))
        return true;
    }
  }
  return false;
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
        continue; // subcluster below threshold: discard, do not call
      }

      string ref = string(chromosome_seqs[chrom] + cl.s, cl.e - cl.s + 1);
      const vector<string> cl_seqs = cl.get_seqs();
      const vector<string> cl_names = cl.get_names();

      // One iteration per POA consensus. Without diploid_poa this runs once with
      // all reads (identical to the previous single-consensus behaviour); with
      // it, a heterozygous cluster yields one iteration per allele, each carrying
      // only its own reads so weight/germline are computed per allele.
      for (const PoaCons &pc : run_poa(cl_seqs)) {
      const int allele_w = (int)pc.read_ids.size();
      if (allele_w < config->min_cluster_weight) {
#pragma omp atomic
        ++n_sub_below_weight;
        continue; // allele below threshold after splitting
      }
      vector<string> allele_names;
      allele_names.reserve(pc.read_ids.size());
      for (int _idx : pc.read_ids)
        if (_idx >= 0 && _idx < (int)cl_names.size())
          allele_names.push_back(cl_names[_idx]);

      vector<SV> _svs;
      const string &consensus = pc.seq;

      // ksw2 stuff - TODO: move to a separate function
      // asm20 scoring (minimap2): match=1 mis=-4 gapo=6,26 gape=2,1.
      // Matches the map-hifi read aligner; avoids asm10's harsh -9 mismatch
      // that fragments single insertions in tandem repeats into spurious pairs.
      int sc_mch = 1, sc_mis = -4, gapo = 6, gape = 2, gapo2 = 26, gape2 = 1;
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
                       allele_w, cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(allele_names);
            _svs.push_back(sv);
            nv++;
          }
          cpos += l;
        } else if (op == 'D') {
          if (l >= config->min_sv_length) {
            SV sv = SV("DEL", cl.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, l),
                       string(chromosome_seqs[chrom] + rpos - 1, 1), allele_w,
                       cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(allele_names);
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
      free(ts);
      free(qs);
      free(ez.cigar); // ksw allocates with km=0 -> standard malloc
      } // end per-consensus loop
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
    // Reciprocal-overlap merge for imprecise (clipped) duplicates of the same
    // event. The length-ratio gate above is too strict for clipped calls, whose
    // breakpoints/lengths are approximate: a real DEL detected from both clip
    // sides can land at slightly different coordinates (length ratio just below
    // min_ratio) and survive as a redundant call. Here we merge same-type
    // imprecise SVs whose reference intervals reciprocally overlap (>= 0.7).
    if (!merged && sv.imprecise && prev.imprecise && sv.type == prev.type &&
        sv.chrom == prev.chrom) {
      long long ov = (long long)min(sv.e, prev.e) - (long long)max(sv.s, prev.s);
      long long un = (long long)max(sv.e, prev.e) - (long long)min(sv.s, prev.s);
      double ro = (un > 0) ? (double)ov / (double)un : 0.0;
      if (ov > 0 && ro >= 0.7) {
        spdlog::debug("[CALLER_CHAIN][RO_MERGE] chrom={} prev_start={} prev_w={} "
                      "new_start={} new_w={} type={} ro={:.2f} keep={}",
                      sv.chrom, prev.s, prev.w, sv.s, sv.w, sv.type, ro,
                      sv.w > prev.w ? "new" : "prev");
        if (sv.w > prev.w)
          prev = sv;
        merged = true;
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
