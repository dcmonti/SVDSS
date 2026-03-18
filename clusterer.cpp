#include "clusterer.hpp"

Clusterer::Clusterer(unordered_map<string, vector<SFS>> *_SFSs, bool _sfs_from_fasta) {
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
      spdlog::debug(
          "[SFS_FILTER][PLACED] chrom={} read={} sfs_qs={} sfs_len={} rs={} re={} qs={} qe={}",
            extsfs.chrom, extsfs.qname, extsfs.qs, extsfs.l, extsfs.rs, extsfs.re, extsfs.qs, extsfs.qe);
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
      spdlog::warn("Non primary alignment.. Is bam file smoothed?");
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
      spdlog::debug(
          "[SFS_FILTER][REV_STRAND] read={} sfs_qs={} sfs_qe={} old_qs={} old_qe={} read len={} sfs len={}",
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
      spdlog::debug(
          "[SFS_FILTER][UNPLACED] read={} chrom={} sfs_qs={} sfs_len={} (both boundaries missing)",
          qname, chrom, sfs.qs, sfs.l);
      ++unplaced;
      continue;
    } else if (refs == -1) {
      uint op = bam_cigar_op(*(cigar + 0));
      uint l = bam_cigar_oplen(*(cigar + 0));
      if (op == BAM_CSOFT_CLIP && config->clipped)
        lclip = make_pair(aln->core.pos, l);
      else {
        spdlog::debug(
            "[SFS_FILTER][START_UNPLACED] read={} chrom={} sfs_qs={} sfs_len={} (no left boundary)",
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
        spdlog::debug(
            "[SFS_FILTER][END_UNPLACED] read={} chrom={} sfs_qs={} sfs_len={} (no right boundary)",
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
      spdlog::warn(
          "[SFS_FILTER][UNKNOWN] read={} chrom={} sfs_qs={} sfs_len={} (failed fallback placement)",
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

  if (lclip.second > 0)
    _p_clips[index].push_back(
        Clip(qname, chrom, lclip.first, lclip.second, true));
  if (rclip.second > 0)
    _p_clips[index].push_back(
        Clip(qname, chrom, rclip.first, rclip.second, false));
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
    spdlog::warn("[CLUSTERING] No extended SFSs available. Skipping clustering.");
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
    spdlog::info(
        "Maximum extended SFS length: {}bp. Using user-defined separation distance: {}bp.",
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
      spdlog::debug(
          "[CLUSTERING][NEW_CHROM] closing interval {}:{}-{} (idx {}-{}; n_sfs={})",
          extended_SFSs[prev_i].chrom, extended_SFSs[prev_i].rs,
          extended_SFSs[i - 1].re, prev_i, i - 1, (i - prev_i));
      prev_chrom = sfs.chrom;
      intervals.push_back(make_pair(prev_i, i - 1));
      prev_i = i;
      prev_e = sfs.re;
      continue;
    } else {
      if (sfs.rs - prev_e > dist) {
        spdlog::debug(
            "[CLUSTERING][DIST_SPLIT] closing interval {}:{}-{} (idx {}-{}; n_sfs={}) gap={} > dist={}",
            extended_SFSs[prev_i].chrom, extended_SFSs[prev_i].rs,
            extended_SFSs[i - 1].re, prev_i, i - 1, (i - prev_i),
            sfs.rs - prev_e, dist);
        intervals.push_back(make_pair(prev_i, i - 1));
        prev_e = sfs.re;
        prev_i = i;
      } else {
        prev_e = max(prev_e, sfs.re);
      }
    }
  }
  intervals.push_back(make_pair(prev_i, extended_SFSs.size() - 1));
  spdlog::debug(
      "[CLUSTERING] Built {} coarse intervals. Last interval {}:{}-{} (idx {}-{}; n_sfs={})",
      intervals.size(), extended_SFSs[prev_i].chrom, extended_SFSs[prev_i].rs,
      extended_SFSs.back().re, prev_i, extended_SFSs.size() - 1,
      (extended_SFSs.size() - prev_i));

  // Cluster SFS inside each interval
  // Use a per-interval vector instead of per-thread maps to avoid
  // thread-dependent clustering results
  vector<map<pair<int, int>, vector<SFS>>> _interval_sfs_clusters(intervals.size());
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < intervals.size(); i++) {
    int j = intervals[i].first;
    int low = extended_SFSs[j].rs;
    int high = extended_SFSs[j].re;
    int last_j = j;
    j++;
    for (; j <= intervals[i].second; j++) {
      const SFS &sfs = extended_SFSs[j];
      if (sfs.rs <= high) {
        low = min(low, sfs.rs);
        high = max(high, sfs.re);
      } else {
        for (int k = last_j; k < j;
             k++) { // CHECKME: < or <=?
                    // NOTE: <= makes the code waaaay slower
          _interval_sfs_clusters[i][make_pair(low, high)].push_back(extended_SFSs[k]);
        }
        low = sfs.rs;
        high = sfs.re;
        last_j = j;
      }
    }
    for (int k = last_j; k <= intervals[i].second;
         k++) { // CHECKME: it was < but in that way
                // we were losing an sfs per cluster
      _interval_sfs_clusters[i][make_pair(low, high)].push_back(extended_SFSs[k]);
    }

    size_t total_sfs = 0;
    for (const auto &cl : _interval_sfs_clusters[i])
      total_sfs += cl.second.size();
    spdlog::debug(
        "[CLUSTERING][INTERVAL_DONE] interval_idx={} chrom={} idx_range=[{}-{}] clusters={} sfs={}",
        i, extended_SFSs[intervals[i].first].chrom, intervals[i].first,
        intervals[i].second, _interval_sfs_clusters[i].size(), total_sfs);
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
  spdlog::debug(
      "[CLUSTERING] Final clusters={} total_assigned_sfs={}.",
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
      spdlog::debug(
          "[CLUSTER_FILTER][LOW_WEIGHT_PRE] chrom={} start={} end={} uniq_reads={} min_required={}",
          cluster.chrom, min_s, max_e, cluster_size, config->min_cluster_weight);
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
        spdlog::debug(
            "[SFS_FILTER][UNEXTENDED_IN_CLUSTER] cluster={}:{}-{} read={} (cannot project full cluster interval on read)",
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
      spdlog::debug(
          "[CLUSTER_FILTER][LOW_WEIGHT_POST] chrom={} start={} end={} subreads={} min_required={}",
          cluster.chrom, min_s, max_e, cluster.size(), config->min_cluster_weight);
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