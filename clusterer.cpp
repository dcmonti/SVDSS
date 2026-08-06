#include "clusterer.hpp"

namespace {

struct SAEntry {
  string chrom;
  uint pos;       // 1-based (SAM spec)
  bool reverse;
  int mapq;
  uint ref_len;
  uint query_start; // offset in SA's own (aligned) orientation
  uint query_len;
};

// Parse a single SA entry "chr,pos,strand,CIGAR,mapQ,NM".
// Returns false if malformed or below min_mapq.
static bool parse_sa_entry(const string &entry, int min_mapq, SAEntry &out) {
  size_t p1 = entry.find(',');
  if (p1 == string::npos) return false;
  size_t p2 = entry.find(',', p1 + 1);
  if (p2 == string::npos) return false;
  size_t p3 = entry.find(',', p2 + 1);
  if (p3 == string::npos) return false;
  size_t p4 = entry.find(',', p3 + 1);
  if (p4 == string::npos) return false;
  size_t p5 = entry.find(',', p4 + 1);
  if (p5 == string::npos) return false;

  out.chrom = entry.substr(0, p1);
  try {
    out.pos = stoul(entry.substr(p1 + 1, p2 - p1 - 1));
    out.mapq = stoi(entry.substr(p4 + 1, p5 - p4 - 1));
  } catch (...) { return false; }
  if (out.mapq < min_mapq) return false;
  out.reverse = (entry[p2 + 1] == '-');

  string sa_cigar = entry.substr(p3 + 1, p4 - p3 - 1);
  out.ref_len = 0;
  out.query_start = 0;
  out.query_len = 0;
  int num = 0;
  bool first_clip_seen = false;
  for (char c : sa_cigar) {
    if (isdigit(c)) {
      num = num * 10 + (c - '0');
    } else {
      if (c == 'S' || c == 'H') {
        if (!first_clip_seen && out.query_len == 0)
          out.query_start = num;
        first_clip_seen = true;
      } else if (c == 'M' || c == '=' || c == 'X') {
        out.ref_len += num;
        out.query_len += num;
        first_clip_seen = true;
      } else if (c == 'I') {
        out.query_len += num;
        first_clip_seen = true;
      } else if (c == 'D' || c == 'N') {
        out.ref_len += num;
        first_clip_seen = true;
      }
      num = 0;
    }
  }
  return true;
}

} // namespace

Clusterer::Clusterer(unordered_map<string, vector<SFS>> *_SFSs,
                     bool _sfs_from_fasta) {
  SFSs = _SFSs;
  sfs_from_fasta = _sfs_from_fasta;
  config = Configuration::getInstance();
}

void Clusterer::run() {
  // 0. Open BAM
  bam_file = hts_open(config->bam.c_str(), "r");
  bam_index = sam_index_load(bam_file, config->bam.c_str());
  bam_header = sam_hdr_read(bam_file);
  bgzf_mt(bam_file->fp.bgzf, 8, 1);

  // 1. Place SFSs on reference genome by extracting subalignments from read
  // alignments and extending them using unique k-mers
  spdlog::info("Placing SFSs on reference genome");
  _p_clips.resize(config->threads);
  _p_extended_sfs.resize(config->threads);
  align_and_extend();
  for (int i = 0; i < config->threads; i++) {
    for (const auto &extsfs : _p_extended_sfs[i]) {
      spdlog::debug("[SFS_FILTER][PLACED] chrom={} read={} sfs_qs={} "
                    "sfs_len={} rs={} re={} qs={} qe={}",
                    extsfs.chrom, extsfs.qname, extsfs.qs, extsfs.l, extsfs.rs,
                    extsfs.re, extsfs.qs, extsfs.qe);
      extended_SFSs.push_back(extsfs);
    }
    clips.insert(clips.begin(), _p_clips[i].begin(), _p_clips[i].end());
  }
  spdlog::info("{}/{}/{} unplaced SFSs. {} erroneus SFSs. {} clipped SFSs.",
               unplaced, s_unplaced, e_unplaced, unknown, clips.size());

  // 2. Cluster SFSs by proximity
  spdlog::info("Clustering {} SFSs..", extended_SFSs.size());
  cluster_by_proximity();
  map<pair<int, int>, Cluster> _ext_clusters;
  for (size_t i = 0; i < _p_sfs_clusters.size(); i++)
    for (const auto &cluster : _p_sfs_clusters[i])
      clusters.push_back(Cluster(cluster.second));

  // 3. Extend SFSs inside each cluster to force them to start/end at same
  // reference position
  spdlog::info("Extending {} clusters..", clusters.size());
  fill_clusters();
  spdlog::info(
      "Filtered {} SFSs. Filtered {} clusters. Filtered {} global clusters.",
      unextended, small_clusters, small_clusters_2);

  // 4. Store clusters to file
  if (config->clusters.compare("") != 0) {
    spdlog::info("Storing clusters to {}", config->clusters);
    store_clusters();
  }

  sam_close(bam_file);
}

/* Extract alignment of SFSs from corresponding read alignment and extend them
 * using unique k-mers in the w-bp flanking regions */
void Clusterer::align_and_extend() {
  // Initializing
  for (int i = 0; i < 2; i++) {
    bam_entries.push_back(vector<vector<bam1_t *>>(config->threads));
    for (int j = 0; j < config->threads; j++)
      for (int k = 0; k < config->batch_size / config->threads; k++)
        bam_entries[i][j].push_back(bam_init1());
  }

  int p = 0;
  int b = 0;
  load_batch(p);
  time_t start_time;
  time_t curr_time;
  time(&start_time);
  bool should_load = true;
  bool should_process = true;
  while (should_process) {
    if (!should_load)
      should_process = false;
#pragma omp parallel for num_threads(config->threads + 1)
    // TODO: avoid additional thread. Make -1 before (same in other classes,
    // with -2)
    for (int i = 0; i < config->threads + 1; i++) {
      int t = omp_get_thread_num();
      if (t == 0) {
        // First thread loads next batch
        if (should_load)
          should_load = load_batch((p + 1) % 2);
      } else
        // Other threads process batch
        process_batch(p, t - 1);
    }
    p += 1;
    p %= 2;
    b += 1;
    time(&curr_time);
    if (curr_time - start_time == 0)
      ++curr_time;
    cerr << "Extended batch " << b << ". Time: " << curr_time - start_time
         << "\r";
  }

  // cleanup
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < config->threads; j++)
      for (int k = 0; k < config->batch_size / config->threads; k++)
        bam_destroy1(bam_entries[i][j][k]);
}

/* Load batch from BAM file and store to input entry p. The logic behind is:
 * fill position i per each thread, then move to position i+1.. */
bool Clusterer::load_batch(int p) {
  int i = 0;
  int nseqs = 0;
  while (sam_read1(bam_file, bam_header,
                   bam_entries[p][nseqs % config->threads][i]) >= 0) {
    bam1_t *aln = bam_entries[p][nseqs % config->threads][i];
    if (aln == nullptr) {
      spdlog::critical("nullptr. Why are we here? Please check");
      exit(1);
    }
    if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
        aln->core.flag & BAM_FSECONDARY) {
      spdlog::debug("Non primary alignment.. Is bam file smoothed?");
      continue;
    }
    if (aln->core.qual < config->min_mapq)
      continue;
    char *qname = bam_get_qname(aln);
    if (SFSs->find(qname) == SFSs->end())
      continue;
    ++nseqs;
    if (nseqs % config->threads == 0)
      ++i;
    if (nseqs == config->batch_size)
      return true;
  }
  // last batch is incomplete since we reached the end of .bam file
  // TODO: can we do like in ping_pong?
  if (nseqs != config->batch_size) {
    for (int j = nseqs % config->threads; j < config->threads; j++)
      for (int _ = i; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
    for (int j = 0; j < nseqs % config->threads; j++)
      for (int _ = i + 1; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
  }
  return false;
}

/* Extend the batch of alignments assigned to thread t */
void Clusterer::process_batch(int p, int t) {
  bam1_t *aln;
  for (size_t b = 0; b < bam_entries[p][t].size(); b++) {
    aln = bam_entries[p][t][b];
    if (aln == nullptr)
      break;
    extend_alignment(aln, t);
  }
}

/* Place SFSs on the alignment and extend them using k-mers */
void Clusterer::extend_alignment(bam1_t *aln, int index) {
  char *qname = bam_get_qname(aln);
  uint32_t *cigar = bam_get_cigar(aln);
  int strand = (aln->core.flag & BAM_FREVERSE) ? -1 : 1;
  vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
  string chrom(bam_header->target_name[aln->core.tid]);
  // TODO: do this check while loading the alignments
  if (chromosome_seqs.find(chrom) == chromosome_seqs.end())
    return;
  // DCM: fix rev compl
  // NOTE: we may have more sfs on a clipped read, but all of them will produce
  // the same clip
  pair<uint, uint> lclip = make_pair(0, 0);
  pair<uint, uint> rclip = make_pair(0, 0);
  int last_pos = 0;
  vector<SFS> local_extended_sfs;
  // If read is on reverse strand AND SFS coordinates come from FASTA
  // (forward-strand coordinates), invert them and re-sort by qs so that the
  // last_pos optimization works correctly. When SFS come from BAM, htslib
  // already provides the sequence in RC orientation, so coordinates match the
  // alignment as-is.
  vector<SFS> &read_sfss = SFSs->at(qname);
  if (strand == -1 && sfs_from_fasta) {
    int read_len = aln->core.l_qseq;
    for (SFS &sfs : read_sfss) {
      int old_qs = sfs.qs;
      int old_qe = sfs.qe;
      sfs.qs = read_len - old_qe;
      sfs.qe = read_len - old_qs;
      spdlog::debug("[SFS_FILTER][REV_STRAND] read={} sfs_qs={} sfs_qe={} "
                    "old_qs={} old_qe={} read len={} sfs len={}",
                    qname, sfs.qs, sfs.qe, old_qs, old_qe, read_len, sfs.l);
    }
    sort(read_sfss.begin(), read_sfss.end(),
         [](const SFS &a, const SFS &b) { return a.qs < b.qs; });
  }
  for (SFS &sfs : read_sfss) {

    int s = sfs.qs;
    int e = sfs.qs + sfs.l - 1;
    int hp_tag = sfs.htag;
    int aln_start = -1;
    int aln_end = -1;
    vector<pair<int, int>> local_alpairs; // subalignment
    // find start and end of SFS in read's alignment
    int refs = -1;
    int refe = -1;
    for (size_t i = last_pos; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      if (q == -1 || r == -1)
        continue;
      else if (q < s) {
        // here we are getting the last "placed" base before the SFS
        // <= seems more correct to me but using < we are more flexible
        last_pos = i;
        refs = r;
        aln_start = i;
      } else if (q > e) {
        // here we are getting the first "placed" base after the SFS
        // >= seems more correct but using > we are more flexible
        refe = r;
        aln_end = i;
        break;
      }
    }
    // Current SFS is aligned from refs to refe. In case of an insertion, the
    // interval will cover the entire insertion

    // We extract the local alignment of the region of interest
    if (refs == -1 && refe == -1) {
      // we couldn't place the first and the last base, so we skip this -
      // otherwise we'll end up considering the entire read
      spdlog::debug("[SFS_FILTER][UNPLACED] read={} chrom={} sfs_qs={} "
                    "sfs_len={} (both boundaries missing)",
                    qname, chrom, sfs.qs, sfs.l);
      ++unplaced;
      continue;
    } else if (refs == -1) {
      uint op = bam_cigar_op(*(cigar + 0));
      uint l = bam_cigar_oplen(*(cigar + 0));
      if (op == BAM_CSOFT_CLIP && config->clipped)
        lclip = make_pair(aln->core.pos, l);
      else {
        spdlog::debug("[SFS_FILTER][START_UNPLACED] read={} chrom={} sfs_qs={} "
                      "sfs_len={} (no left boundary)",
                      qname, chrom, sfs.qs, sfs.l);
        ++s_unplaced;
      }
      continue; // in any case, we skip this SFS
    } else if (refe == -1) {
      uint op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
      uint l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
      if (op == BAM_CSOFT_CLIP && config->clipped)
        rclip = make_pair(bam_endpos(aln), l);
      else {
        spdlog::debug("[SFS_FILTER][END_UNPLACED] read={} chrom={} sfs_qs={} "
                      "sfs_len={} (no right boundary)",
                      qname, chrom, sfs.qs, sfs.l);
        ++e_unplaced;
      }
      continue; // in any case, we skip this SFS
    } else {
      // we placed the first and last base, so we extract the subalignment
      int last_r = refs - 1;
      for (int i = aln_start; i <= aln_end; i++) {
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (r == -1) {
          if (refs <= last_r && last_r <= refe)
            local_alpairs.push_back(make_pair(q, r));
        } else {
          last_r = r;
          if (refs <= r && r <= refe)
            local_alpairs.push_back(make_pair(q, r));
        }
        // We break when we found a placed base at or after the reference end
        if (q != -1 && r != -1 && r >= refe)
          break;
      }
    }

    // SFS has been placed and local_alpairs contains the subalignment

    // extract the config->flank pairs preceding the region of interest
    vector<pair<int, int>> pre_alpairs;
    uint n = 0;
    for (int i = aln_start - 1; i >= 0; --i) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      pre_alpairs.push_back(make_pair(q, r));
      ++n;
      if (n == config->flank)
        break;
    }
    reverse(pre_alpairs.begin(), pre_alpairs.end());
    // extract the config->flank pairs following the region of interest
    vector<pair<int, int>> post_alpairs;
    n = 0;
    for (uint i = aln_end + 1; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      post_alpairs.push_back(make_pair(q, r));
      ++n;
      if (n == config->flank)
        break;
    }

    // get the unique kmer in the upstream and downstream config->flank bp
    // regions
    pair<int, int> prekmer =
        get_unique_kmers(pre_alpairs, config->ksize, true,
                         chrom); // true for first kmer found (shorter cluster)
    pair<int, int> postkmer =
        get_unique_kmers(post_alpairs, config->ksize, false,
                         chrom); // false for first kmer found (shorter cluster)

    // if we couldn't place a kmer, we just get the entire region
    if (prekmer.first == -1 || prekmer.second == -1) {
      prekmer.first = local_alpairs.front().first;
      prekmer.second = local_alpairs.front().second;
    }
    if (postkmer.first == -1 || postkmer.second == -1) {
      postkmer.first = local_alpairs.back().first;
      postkmer.second = local_alpairs.back().second;
    }
    // if also the entire region is not correctly placed, then we skip it
    // NOTE: I think we can solve this by increasing config->flank.
    if (prekmer.first == -1 || prekmer.second == -1 || postkmer.first == -1 ||
        postkmer.second == -1) {
      spdlog::warn("[SFS_FILTER][UNKNOWN] read={} chrom={} sfs_qs={} "
                   "sfs_len={} (failed fallback placement)",
                   qname, chrom, sfs.qs, sfs.l);
      ++unknown;
      continue;
    }
    // FIXME: understand why this is happening (chr16 on full giab genome)
    if ((uint)prekmer.second > postkmer.second + config->ksize) {
      spdlog::warn("Error on {}. SFS starting at {} (length {})", qname, sfs.qs,
                   sfs.l);
    } else {
      SFS ext_sfs(string(chrom), string(qname), prekmer.second,
                  postkmer.second + config->ksize, prekmer.first,
                  postkmer.first + config->ksize, hp_tag);
      ext_sfs.add_orig_interval(refs, refe);
      local_extended_sfs.push_back(ext_sfs);
    }
  }
  // When two SFSs are close but not overlapping, we may end up with two
  // overlapping extended SFSs. We need to merge SFSs two by two until we don't
  // need to merge anything
  vector<SFS> merged_extended_sfs;
  for (size_t i = 0; i < local_extended_sfs.size(); ++i) {
    size_t j;
    for (j = 0; j < merged_extended_sfs.size(); ++j) {
      if ((local_extended_sfs.at(i).rs <= merged_extended_sfs.at(j).rs &&
           merged_extended_sfs.at(j).rs <= local_extended_sfs.at(i).re) ||
          (merged_extended_sfs.at(j).rs <= local_extended_sfs.at(i).rs &&
           local_extended_sfs.at(i).rs <= merged_extended_sfs.at(j).re))
        break;
    }
    if (j < merged_extended_sfs.size()) {
      merged_extended_sfs[j].rs =
          min(merged_extended_sfs.at(j).rs, local_extended_sfs.at(i).rs);
      merged_extended_sfs[j].re =
          max(merged_extended_sfs.at(j).re, local_extended_sfs.at(i).re);
      merged_extended_sfs[j].qs =
          min(merged_extended_sfs.at(j).qs, local_extended_sfs.at(i).qs);
      merged_extended_sfs[j].qe =
          max(merged_extended_sfs.at(j).qe, local_extended_sfs.at(i).qe);
      merged_extended_sfs[j].orig_intervals.insert(
          merged_extended_sfs[j].orig_intervals.end(),
          local_extended_sfs.at(i).orig_intervals.begin(),
          local_extended_sfs.at(i).orig_intervals.end());
    } else {
      merged_extended_sfs.push_back(local_extended_sfs.at(i));
    }
  }
  for (const auto &mes : merged_extended_sfs)
    _p_extended_sfs[index].push_back(mes);

  // Parse all SA entries from the SA tag (one read may have multiple
  // supplementary alignments, e.g. complex rearrangements).
  bool primary_reverse = (aln->core.flag & BAM_FREVERSE) != 0;
  vector<SAEntry> sa_entries;
  uint8_t *sa_tag = bam_aux_get(aln, "SA");
  if (sa_tag != NULL) {
    string sa_str(bam_aux2Z(sa_tag));
    size_t s = 0;
    while (s < sa_str.size()) {
      size_t e = sa_str.find(';', s);
      if (e == string::npos) break;
      if (e > s) {
        SAEntry sae;
        if (parse_sa_entry(sa_str.substr(s, e - s), config->min_mapq, sae))
          sa_entries.push_back(sae);
      }
      s = e + 1;
    }
  }

  // Compute primary's aligned query span from CIGAR (independent of SFS).
  uint read_len = aln->core.l_qseq;
  uint primary_q_start = 0;
  uint primary_q_end = read_len;
  if (aln->core.n_cigar > 0) {
    uint first_op = bam_cigar_op(*(cigar + 0));
    uint first_l = bam_cigar_oplen(*(cigar + 0));
    if (first_op == BAM_CSOFT_CLIP || first_op == BAM_CHARD_CLIP)
      primary_q_start = first_l;
    uint last_op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
    uint last_l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
    if (last_op == BAM_CSOFT_CLIP || last_op == BAM_CHARD_CLIP)
      primary_q_end = (read_len > last_l) ? (read_len - last_l) : 0;
  }

  // Convert an SA's query span into the primary's read orientation, so we can
  // compare it against primary_q_start / primary_q_end directly. SA CIGAR is
  // expressed in SA's own aligned orientation; if SA strand differs from primary,
  // its query coordinates need to be flipped onto the primary frame.
  auto sa_qspan_in_primary = [&](const SAEntry &sa) -> pair<uint, uint> {
    if (sa.reverse == primary_reverse)
      return {sa.query_start, sa.query_start + sa.query_len};
    uint sa_end = sa.query_start + sa.query_len;
    uint flipped_end = (read_len > sa.query_start) ? (read_len - sa.query_start) : 0;
    uint flipped_start = (read_len > sa_end) ? (read_len - sa_end) : 0;
    return {flipped_start, flipped_end};
  };

  // Fill in this read's junction geometry: dq (inner unaligned query bases) and,
  // when dq > 0, the inserted bases themselves.
  //
  // dq must be computed from query coordinates in the PRIMARY's frame, which is
  // what sa_qspan_in_primary yields — the Clip constructor's fallback uses the
  // raw SA CIGAR offsets and is therefore only valid for a same-strand split.
  //
  //   left  clip: the clip occupies primary-frame query [0, l);
  //               the SA covers [span.first, span.second) -> dq = l - span.second
  //   right clip: the clip occupies [primary_q_end, primary_q_end + l);
  //               the SA starts at span.first          -> dq = span.first - primary_q_end
  //
  // Sign is meaningful and must survive: dq > 0 is novel sequence between the
  // segments (SVINSLEN/SVINSSEQ), dq < 0 is breakpoint microhomology
  // (HOMLEN/HOMSEQ, recoverable from the reference so no bases are stored here).
  auto set_junction = [&](Clip &c, const SAEntry &sa, uint clip_len,
                          bool is_left) {
    auto span = sa_qspan_in_primary(sa);
    long long dq = is_left ? (long long)clip_len - (long long)span.second
                           : (long long)span.first - (long long)primary_q_end;
    c.dq = (int)dq;
    c.dqs.assign(1, c.dq);
    if (dq <= 0)
      return; // microhomology or blunt: nothing novel to store
    // Novel bases, in primary-frame query coordinates. bam_get_seq is stored in
    // the primary's orientation too, so the two frames agree.
    long long qs = is_left ? (long long)span.second : (long long)primary_q_end;
    if (qs < 0 || qs + dq > (long long)read_len)
      return;
    const uint8_t *seq = bam_get_seq(aln);
    string ins;
    ins.reserve((size_t)dq);
    for (long long k = 0; k < dq; k++)
      ins.push_back(seq_nt16_str[bam_seqi(seq, qs + k)]);
    c.ins_seq = std::move(ins);
  };

  // Pick the SA whose query span is most adjacent to the clip on the read.
  // For a left clip we match SA's end against primary's aligned start; for a
  // right clip we match SA's start against primary's aligned end. Tiebreakers
  // (in order): higher mapq, longer SA reference span.
  auto pick_best_sa = [&](uint target, bool match_end) -> int {
    int best = -1;
    long long best_dist = 0;
    for (size_t k = 0; k < sa_entries.size(); k++) {
      auto span = sa_qspan_in_primary(sa_entries[k]);
      long long d = match_end
          ? abs((long long)span.second - (long long)target)
          : abs((long long)span.first - (long long)target);
      if (best < 0 || d < best_dist) {
        best = (int)k;
        best_dist = d;
      } else if (d == best_dist) {
        const SAEntry &b = sa_entries[best];
        const SAEntry &c = sa_entries[k];
        if (c.mapq > b.mapq ||
            (c.mapq == b.mapq && c.ref_len > b.ref_len))
          best = (int)k;
      }
    }
    return best;
  };

  if (lclip.second > 0) {
    int best = sa_entries.empty() ? -1 : pick_best_sa(primary_q_start, true);
    if (best >= 0) {
      const SAEntry &sa = sa_entries[best];
      Clip c(qname, chrom, lclip.first, lclip.second, true, primary_reverse,
             sa.reverse, sa.chrom, sa.pos, sa.ref_len, sa.query_start,
             sa.query_len);
      set_junction(c, sa, lclip.second, true);
      _p_clips[index].push_back(std::move(c));
    } else {
      _p_clips[index].push_back(
          Clip(qname, chrom, lclip.first, lclip.second, true));
    }
  }
  if (rclip.second > 0) {
    int best = sa_entries.empty() ? -1 : pick_best_sa(primary_q_end, false);
    if (best >= 0) {
      const SAEntry &sa = sa_entries[best];
      Clip c(qname, chrom, rclip.first, rclip.second, false, primary_reverse,
             sa.reverse, sa.chrom, sa.pos, sa.ref_len, sa.query_start,
             sa.query_len);
      set_junction(c, sa, rclip.second, false);
      _p_clips[index].push_back(std::move(c));
    } else {
      _p_clips[index].push_back(
          Clip(qname, chrom, rclip.first, rclip.second, false));
    }
  }
}

/* Get first/last kmer entirely mapped and with a single occurrence in a
 * subalignment of interest. Return the kmer as the initial pair of positions in
 * the alignment */
pair<int, int>
Clusterer::get_unique_kmers(const vector<pair<int, int>> &alpairs, const uint k,
                            const bool from_end, string chrom) {
  if (alpairs.size() < k)
    return make_pair(-1, -1);

  // Do kmer counting in region
  map<string, int> kmers;
  string kmer_seq;
  size_t i = 0;
  while (i < alpairs.size() - k + 1) {
    bool skip = false;
    for (size_t j = i; j < i + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        // we want clean kmers only - ie placed kmers, no insertions or
        // deletions
        skip = true;
        i = j + 1; // jump to next clean/placed position
        break;
      }
    }
    if (skip)
      continue;
    string kmer(chromosome_seqs[chrom] + alpairs[i].second, k);
    ++kmers[kmer];
    ++i;
  }
  // Get first/last kmer with single occurrence
  pair<int, int> last_kmer = make_pair(-1, -1);
  i = 0;
  while (i < alpairs.size() - k + 1) {
    int offset = i;
    if (from_end)
      offset = alpairs.size() - k - i;
    assert(offset >= 0);
    bool skip = false;
    for (size_t j = offset; j < offset + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        skip = true;
        i += (j - offset); // jump to next possible start position
        break;
      }
    }
    if (skip) {
      ++i;
      continue;
    }
    last_kmer = alpairs[offset];
    string kmer(chromosome_seqs[chrom] + alpairs[offset].second, k);
    if (kmers[kmer] == 1)
      break;
    ++i;
  }
  return last_kmer;
}

void Clusterer::cluster_by_proximity() {
  if (extended_SFSs.empty()) {
    spdlog::warn(
        "[CLUSTERING] No extended SFSs available. Skipping clustering.");
    _p_sfs_clusters.clear();
    return;
  }

  sort(extended_SFSs.begin(), extended_SFSs.end());
  auto r = max_element(extended_SFSs.begin(), extended_SFSs.end(),
                       [](const SFS &lhs, const SFS &rhs) {
                         return lhs.re - lhs.rs < rhs.re - rhs.rs;
                       });
  int dist;
  if (config->max_cluster_dist > 0) {
    dist = config->max_cluster_dist;
    spdlog::info("Maximum extended SFS length: {}bp. Using user-defined "
                 "separation distance: {}bp.",
                 r->re - r->rs, dist);
  } else {
    dist = (r->re - r->rs) * 1.1;
    spdlog::info(
        "Maximum extended SFS length: {}bp. Using separation distance: {}bp.",
        r->re - r->rs, dist);
  }
  spdlog::debug(
      "[CLUSTERING] Start proximity clustering on {} SFSs with distance={}.",
      extended_SFSs.size(), dist);

  // Cluster SFSs inside dist-bp windows
  int prev_i = 0;
  int prev_e = extended_SFSs[0].re;
  string prev_chrom = extended_SFSs[0].chrom;
  vector<pair<int, int>> intervals;
  for (size_t i = 1; i < extended_SFSs.size(); i++) {
    const auto &sfs = extended_SFSs[i];
    // new chromosome
    if (sfs.chrom != prev_chrom) {
      spdlog::debug("[CLUSTERING][NEW_CHROM] closing interval {}:{}-{} (idx "
                    "{}-{}; n_sfs={})",
                    extended_SFSs[prev_i].chrom, extended_SFSs[prev_i].rs,
                    extended_SFSs[i - 1].re, prev_i, i - 1, (i - prev_i));
      prev_chrom = sfs.chrom;
      intervals.push_back(make_pair(prev_i, i - 1));
      prev_i = i;
      prev_e = sfs.re;
      continue;
    } else {
      if (sfs.rs - prev_e > dist) {
        spdlog::debug("[CLUSTERING][DIST_SPLIT] closing interval {}:{}-{} (idx "
                      "{}-{}; n_sfs={}) gap={} > dist={}",
                      extended_SFSs[prev_i].chrom, extended_SFSs[prev_i].rs,
                      extended_SFSs[i - 1].re, prev_i, i - 1, (i - prev_i),
                      sfs.rs - prev_e, dist);
        intervals.push_back(make_pair(prev_i, i - 1));
        prev_e = sfs.re;
        prev_i = i;
      }
      // else {
      //   prev_e = max(prev_e, sfs.re);
      // }
    }
  }
  intervals.push_back(make_pair(prev_i, extended_SFSs.size() - 1));
  spdlog::debug("[CLUSTERING] Built {} coarse intervals. Last interval "
                "{}:{}-{} (idx {}-{}; n_sfs={})",
                intervals.size(), extended_SFSs[prev_i].chrom,
                extended_SFSs[prev_i].rs, extended_SFSs.back().re, prev_i,
                extended_SFSs.size() - 1, (extended_SFSs.size() - prev_i));

  // Cluster SFS inside each interval
  // Use a per-interval vector instead of per-thread maps to avoid
  // thread-dependent clustering results
  vector<map<pair<int, int>, vector<SFS>>> _interval_sfs_clusters(
      intervals.size());
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < intervals.size(); i++) {
    int start_j = intervals[i].first;
    int end_j = intervals[i].second;
    int n_sfs = end_j - start_j + 1;

    // Base case: we have 0 or 1 SFS, cluster is already done
    if (n_sfs == 0) {
      continue;
    } else if (n_sfs == 1) {
      const SFS &sfs = extended_SFSs[start_j];
      _interval_sfs_clusters[i][make_pair(sfs.rs, sfs.re)].push_back(sfs);
      continue;
    }

    // Start with hierarchical agglomerative clustering:
    // 1. Build a similarity distance for each sfs pair A,B in this cluster,
    // using the formula:
    //    1 - RO(A, B), with RO(A, B) = |[A.rs, A.re] ^ [B.rs, B.re]| /
    //    max(A.re-A.rs, B.re-B.rs)
    //    TODO: I'm using the original SFS coords pre extension, is it
    //    necessary?
    // 2. Each SFS is a cluster, iterate and merge the two most similar clusters
    // 3. STOP criteria: when we have a big decrease in similarity.
    //    we are using Mojena's rule: iter[k] - iter[k+1] > mu + c . sigma (mean
    //    and std of the different iterations)

    // init
    struct LocalCluster {
      vector<int> sfs_indices;
      bool active = true;
    };
    vector<LocalCluster> local_clusters(n_sfs);
    for (int k = 0; k < n_sfs; ++k) {
      local_clusters[k].sfs_indices.push_back(k);
    }

    // similarity distance
    vector<vector<float>> dist(n_sfs, vector<float>(n_sfs, 0.0f));
    for (int u = 0; u < n_sfs; ++u) {
      for (int v = u + 1; v < n_sfs; ++v) {
        const SFS &a = extended_SFSs[start_j + u];
        const SFS &b = extended_SFSs[start_j + v];

        int a_orig_rs = numeric_limits<int>::max();
        int a_orig_re = 0;
        for (const auto &interval : a.orig_intervals) {
          a_orig_rs = min(a_orig_rs, interval.first);
          a_orig_re = max(a_orig_re, interval.second);
        }

        int b_orig_rs = numeric_limits<int>::max();
        int b_orig_re = 0;
        for (const auto &interval : b.orig_intervals) {
          b_orig_rs = min(b_orig_rs, interval.first);
          b_orig_re = max(b_orig_re, interval.second);
        }

        int overlap =
            max(0, min(a_orig_re, b_orig_re) - max(a_orig_rs, b_orig_rs));
        float len_a = (float)(a_orig_re - a_orig_rs);
        float len_b = (float)(b_orig_re - b_orig_rs);

        float max_len = max(len_a, len_b);
        float ro = (max_len > 0) ? ((float)overlap / max_len) : 0.0f;

        dist[u][v] = dist[v][u] = 1.0f - ro;
      }
    }

    // clustering with trace
    struct MergeState {
      float merge_distance;
      vector<vector<int>> clusters;
    };
    vector<MergeState> states;

    vector<vector<int>> initial_state;
    for (int k = 0; k < n_sfs; ++k)
      initial_state.push_back({k});
    states.push_back({0.0f, initial_state});

    int active_count = n_sfs;
    while (active_count > 1) {
      float min_d = 2.0f; //
      int best_u = -1, best_v = -1;

      for (int u = 0; u < n_sfs; ++u) {
        if (!local_clusters[u].active)
          continue;
        for (int v = u + 1; v < n_sfs; ++v) {
          if (!local_clusters[v].active)
            continue;
          if (dist[u][v] < min_d) {
            min_d = dist[u][v];
            best_u = u;
            best_v = v;
          }
        }
      }

      if (best_u == -1 || min_d > 0.99f)
        break;

      local_clusters[best_u].sfs_indices.insert(
          local_clusters[best_u].sfs_indices.end(),
          local_clusters[best_v].sfs_indices.begin(),
          local_clusters[best_v].sfs_indices.end());
      local_clusters[best_v].active = false;

      for (int k = 0; k < n_sfs; ++k) {
        if (k != best_u && local_clusters[k].active) {
          float new_dist = max(dist[best_u][k], dist[best_v][k]);
          dist[best_u][k] = dist[k][best_u] = new_dist;
        }
      }
      active_count--;

      vector<vector<int>> current_clusters;
      for (int k = 0; k < n_sfs; ++k) {
        if (local_clusters[k].active) {
          current_clusters.push_back(local_clusters[k].sfs_indices);
        }
      }
      states.push_back({min_d, current_clusters});
    }

    // find max jump in similarity
    int best_state_idx = states.size() - 1; 
    float max_jump = -1.0f;

    for (size_t s = 1; s < states.size(); ++s) {
      float jump = states[s].merge_distance - states[s - 1].merge_distance;
      if (jump > max_jump && states[s].merge_distance > 0.4f) {
        max_jump = jump;
        best_state_idx =
            s - 1; 
      }
    }

    for (const auto &cl_indices : states[best_state_idx].clusters) {
      int low = extended_SFSs[start_j + cl_indices[0]].rs;
      int high = extended_SFSs[start_j + cl_indices[0]].re;

      for (int idx : cl_indices) {
        low = min(low, extended_SFSs[start_j + idx].rs);
        high = max(high, extended_SFSs[start_j + idx].re);
      }

      for (int idx : cl_indices) {
        _interval_sfs_clusters[i][make_pair(low, high)].push_back(
            extended_SFSs[start_j + idx]);
      }
    }
    size_t total_sfs = 0;
    for (const auto &cl : _interval_sfs_clusters[i])
      total_sfs += cl.second.size();
    spdlog::debug("[CLUSTERING][INTERVAL_DONE_HAC] interval_idx={} chrom={} "
                  "idx_range=[{}-{}] clusters={} sfs={}",
                  i, extended_SFSs[intervals[i].first].chrom,
                  intervals[i].first, intervals[i].second,
                  _interval_sfs_clusters[i].size(), total_sfs);
  }
  // Flatten into _p_sfs_clusters for compatibility
  _p_sfs_clusters.resize(1);
  for (size_t i = 0; i < intervals.size(); i++)
    for (const auto &cluster : _interval_sfs_clusters[i])
      _p_sfs_clusters[0][cluster.first] = cluster.second;

  size_t total_clusters = _p_sfs_clusters[0].size();
  size_t total_sfs = 0;
  for (const auto &cluster : _p_sfs_clusters[0])
    total_sfs += cluster.second.size();
  spdlog::info("[CLUSTERING] Final clusters={} total_assigned_sfs={}.",
                total_clusters, total_sfs);
}

// /* Assign coverage and read (sub)sequence to each cluster  */
void Clusterer::fill_clusters() {
  // Allocate
  char *seq[config->threads];
  uint32_t len[config->threads];
  bam1_t *_p_aln[config->threads];
  samFile *_p_bam_file[config->threads];
  hts_idx_t *_p_bam_index[config->threads];
  bam_hdr_t *_p_bam_header[config->threads];
  for (int i = 0; i < config->threads; i++) {
    len[i] = 0;
    _p_aln[i] = bam_init1();
    _p_bam_file[i] = hts_open(config->bam.c_str(), "r");
    _p_bam_index[i] = sam_index_load(_p_bam_file[i], config->bam.c_str());
    _p_bam_header[i] = sam_hdr_read(_p_bam_file[i]);
    bgzf_mt(_p_bam_file[i]->fp.bgzf, 8, 1);
  }

#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < clusters.size(); i++) {
    int t = omp_get_thread_num();
    Cluster &cluster = clusters[i];

    // Force all extended SFSs to start and end at the same position. Build a
    // "global" cluster
    set<string> reads;
    int min_s = numeric_limits<int>::max();
    int max_e = 0;
    for (const SFS &sfs : cluster.SFSs) {
      min_s = min(min_s, sfs.rs);
      max_e = max(max_e, sfs.re);
      reads.insert(sfs.qname);
    }

    size_t cluster_size = reads.size();
    if (cluster_size < config->min_cluster_weight) {
      spdlog::debug("[CLUSTER_FILTER][LOW_WEIGHT_PRE] chrom={} start={} end={} "
                    "uniq_reads={} min_required={}",
                    cluster.chrom, min_s, max_e, cluster_size,
                    config->min_cluster_weight);
      ++small_clusters;
      continue;
    }

    cluster.set_coordinates(min_s, max_e);

    // Iterate over alignments falling in the cluster region to: (i) get total
    // number of reads and (ii) get SFS sequence, one per read
    vector<int> coverages(3, 0);
    vector<tuple<int, int>> locus_reads;

    string region =
        cluster.chrom + ":" + to_string(min_s) + "-" + to_string(max_e);
    hts_itr_t *itr =
        sam_itr_querys(_p_bam_index[t], _p_bam_header[t], region.c_str());
    while (sam_itr_next(_p_bam_file[t], itr, _p_aln[t]) > 0) {
      bam1_t *aln = _p_aln[t];
      if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
          aln->core.flag & BAM_FSECONDARY)
        continue;
      if (aln->core.qual < config->min_mapq)
        continue;
      int hp_t = bam_aux_get(aln, "HP") != NULL
                     ? bam_aux2i(bam_aux_get(aln, "HP"))
                     : 0;
      ++coverages[hp_t]; // FIXME: this cov takes into account also reads
                         // starting or ending inside the cluster (maybe we
                         // should skip those?)
      locus_reads.push_back(make_tuple(0, hp_t == 0 ? 3 : hp_t));
      char *qname = bam_get_qname(aln);
      if (reads.find(qname) == reads.end())
        // we do not have a SFS on this read at this locus
        continue;
      get<0>(locus_reads.back()) = 1;

      // If we have a SFS on this read, load the read sequence
      uint32_t l = aln->core.l_qseq;
      if (l >= len[t]) {
        if (len[t] != 0)
          free(seq[t]);
        len[t] = l;
        seq[t] = (char *)malloc(l + 1);
      }
      uint8_t *q = bam_get_seq(aln);
      for (uint i = 0; i < l; i++)
        seq[t][i] = seq_nt16_str[bam_seqi(q, i)];
      seq[t][l] = '\0';

      // Extract "global" sequence by getting start/end position on read
      // sequence aligning to start/end of cluster
      vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
      int qs = -1, qe = -1;
      for (int i = alpairs.size() - 1; i >= 0; --i) {
        // finding starting position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1)
          continue;
        if (r <= min_s) {
          qs = q;
          break;
        }
      }
      for (uint i = 0; i < alpairs.size(); ++i) {
        // finding ending position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1)
          continue;
        if (r >= max_e) {
          qe = q;
          break;
        }
      }
      if (qs == -1 || qe == -1) {
        // reads starts or ends inside the cluster
        // TODO: get only remaining prefix/suffix? but this may make POA and
        // realignment harder
        spdlog::debug("[SFS_FILTER][UNEXTENDED_IN_CLUSTER] cluster={}:{}-{} "
                      "read={} (cannot project full cluster interval on read)",
                      cluster.chrom, min_s, max_e, qname);
        ++unextended;
      } else {
        string _seq(seq[t], qs, qe - qs + 1);
        cluster.add_subread(qname, _seq, hp_t);
      }
    }
    if (cluster.size() >= config->min_cluster_weight) {
      cluster.set_cov(coverages);
      cluster.set_reads(locus_reads);
    } else {
      spdlog::debug("[CLUSTER_FILTER][LOW_WEIGHT_POST] chrom={} start={} "
                    "end={} subreads={} min_required={}",
                    cluster.chrom, min_s, max_e, cluster.size(),
                    config->min_cluster_weight);
      ++small_clusters_2;
    }
  }

  // clean
  for (int i = 0; i < config->threads; i++) {
    if (len[i] > 0)
      free(seq[i]);
    bam_destroy1(_p_aln[i]);
    sam_close(_p_bam_file[i]);
  }
}

/* Store clusters to file */
void Clusterer::store_clusters() {
  ofstream clofile;
  clofile.open(config->clusters);
  for (const auto &cluster : clusters) {
    // CHECKME: is ending position (cluster.first.second) inclusive? I'm
    // assuming yes
    clofile << cluster.chrom << ":" << cluster.s + 1 << "-" << cluster.e + 1
            << "\t" << cluster.size();
    for (size_t i = 0; i < cluster.size(); ++i)
      clofile << "\t" << cluster.get_name(i) << ":" << cluster.get_seq(i);
    clofile << endl;
  }
  clofile.close();
}