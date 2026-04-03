#include "clipper.hpp"

Clipper::Clipper(const vector<Clip> &_clips) { clips = _clips; }

vector<Clip> Clipper::remove_duplicates(const vector<Clip> &clips) {
  vector<Clip> unique_clips;
  unordered_map<string, int> qnames;
  for (const Clip &clip : clips) {
    if (qnames.find(clip.name) == qnames.end()) {
      qnames[clip.name] = 0;
      unique_clips.push_back(clip);
    }
  }
  return unique_clips;
}

vector<Clip> Clipper::combine(const vector<Clip> &clips) {
  int threads = 4;
  vector<vector<Clip>> _p_combined_clips;
  _p_combined_clips.resize(threads);
  // we first cluster by breakpoints
  unordered_map<string, unordered_map<uint, vector<Clip>>> clips_dict;
  for (const Clip &c : clips) {
    clips_dict[c.chrom][c.p].push_back(c);
  }
// we then merge
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (uint i = 0; i < chromosomes.size(); i++) {
    int t = i % threads;
    const string &chrom = chromosomes[i];
    for (auto it = clips_dict[chrom].begin(); it != clips_dict[chrom].end();
         ++it) {
      uint max_l = 0;
      bool has_sa = false;
      string sa_chrom = "";
      uint sa_pos = 0;
      for (const Clip &c : it->second) {
        if (c.l > max_l) {
          max_l = c.l;
        }
        if (c.sa_has_info) {
          has_sa = true;
          sa_chrom = c.sa_chrom;
          sa_pos = c.sa_pos;
        }
      }
      Clip clip = Clip("", chrom, it->first, max_l, it->second.front().starting,
                       it->second.size());
      if (has_sa) {
        clip.sa_chrom = sa_chrom;
        clip.sa_pos = sa_pos;
        clip.sa_has_info = true;
      }
      _p_combined_clips[t].push_back(clip);
    }
  }
  vector<Clip> combined_clips;
  for (int i = 0; i < threads; i++) {
    combined_clips.insert(combined_clips.begin(), _p_combined_clips[i].begin(),
                          _p_combined_clips[i].end());
  }
  return combined_clips;
}

vector<Clip> Clipper::filter_lowcovered(const vector<Clip> &clips,
                                        const uint w) {
  vector<Clip> filtered_clips;
  for (const Clip &c : clips) {
    if (c.w >= w) {
      filtered_clips.push_back(c);
    }
  }
  return filtered_clips;
}

// Cluster clips by proximity
// TODO: this might be too slow
vector<Clip> Clipper::cluster(const vector<Clip> &clips, uint r) {
  vector<Clip> clusters;
  map<uint, Clip> clusters_by_pos;
  for (const Clip &c : clips) {
    bool found = false;
    for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
         it != clusters_by_pos.end(); ++it) {
      if (it->first - r <= c.p && c.p <= it->first + r) {
        found = true;
        it->second.l = max(it->second.l, c.l);
        it->second.w += c.w;
        if (c.sa_has_info) {
          it->second.sa_chrom = c.sa_chrom;
          it->second.sa_pos = c.sa_pos;
          it->second.sa_has_info = true;
        }
      }
    }
    if (!found) {
      clusters_by_pos[c.p] = c;
    }
  }

  for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
       it != clusters_by_pos.end(); ++it) {
    clusters.push_back(it->second);
  }
  return clusters;
}

vector<Clip> Clipper::filter_tooclose_clips(const vector<Clip> &clips,
                                            interval_tree_t<int> &vartree) {
  vector<Clip> fclips;
  for (const Clip &c : clips) {
    if (vartree.overlap_find({(int)c.p, (int)c.p + 1}) == end(vartree)) {
      fclips.push_back(c);
    }
  }
  return fclips;
}

// find smallest right that is larger than query
int binary_search(const vector<Clip> &clips, uint begin, uint end,
                  const Clip &query) {
  if (begin > end || begin >= clips.size()) {
    return -1;
  }
  uint m = (begin + end) / 2;
  if (m >= clips.size()) {
    return -1;
  }
  if (clips[m].p == query.p) {
    if (m + 1 < clips.size()) {
      return m + 1;
    } else {
      return m;
    }
  } else if (clips[m].p > query.p) {
    if (m > 0 && clips[m - 1].p < query.p) {
      return m;
    }
    if (m == 0) {
      return 0; // query is before all clips, return first element
    }
    return binary_search(clips, begin, m - 1, query);
  } else {
    return binary_search(clips, m + 1, end, query);
  }
}

void Clipper::call(int threads, interval_tree_t<int> &vartree) {
  vector<Clip> rclips;
  vector<Clip> lclips;
  for (const Clip &clip : clips) {
    if (clip.starting) {
      lclips.push_back(clip);
    } else {
      rclips.push_back(clip);
    }
  }
  spdlog::info("Clipped SFS: {} left clips, {} right clips.", lclips.size(),
               rclips.size());
#pragma omp parallel for num_threads(2) schedule(static, 1)
  for (int i = 0; i < 2; i++) {
    if (i == 0) {
      rclips = remove_duplicates(rclips);
      spdlog::info("[CLIP_FILTER][RIGHT] after remove_duplicates: {}", rclips.size());
      rclips = combine(rclips);
      spdlog::info("[CLIP_FILTER][RIGHT] after combine: {}", rclips.size());
      rclips = filter_lowcovered(rclips, 1);
      spdlog::info("[CLIP_FILTER][RIGHT] after filter_lowcovered(1): {}", rclips.size());
      rclips = filter_tooclose_clips(rclips, vartree);
      spdlog::info("[CLIP_FILTER][RIGHT] after filter_tooclose_clips: {}", rclips.size());
      rclips = cluster(rclips, 1000);
      spdlog::info("[CLIP_FILTER][RIGHT] after cluster(1000): {}", rclips.size());
      sort(rclips.begin(), rclips.end());
    } else {
      lclips = remove_duplicates(lclips);
      spdlog::info("[CLIP_FILTER][LEFT] after remove_duplicates: {}", lclips.size());
      lclips = combine(lclips);
      spdlog::info("[CLIP_FILTER][LEFT] after combine: {}", lclips.size());
      lclips = filter_lowcovered(lclips, 1);
      spdlog::info("[CLIP_FILTER][LEFT] after filter_lowcovered(1): {}", lclips.size());
      lclips = filter_tooclose_clips(lclips, vartree);
      spdlog::info("[CLIP_FILTER][LEFT] after filter_tooclose_clips: {}", lclips.size());
      lclips = cluster(lclips, 1000);
      spdlog::info("[CLIP_FILTER][LEFT] after cluster(1000): {}", lclips.size());
      sort(lclips.begin(), lclips.end());
    }
  }
  spdlog::info("After filtering: {} left clips, {} right clips.", lclips.size(),
               rclips.size());
  _p_svs.resize(threads);
  if (lclips.empty() || rclips.empty()) {
    return;
  }
  // Predicting insertions
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (uint i = 0; i < lclips.size(); i++) {
    const Clip &lc = lclips[i];
    int t = omp_get_thread_num();
    string chrom = lc.chrom;
    
    bool sa_used = false;
    if (lc.sa_has_info && lc.sa_chrom == chrom) {
       uint s = min(lc.p, lc.sa_pos);
       uint e = max(lc.p, lc.sa_pos);
       uint l = e - s;
       if (l < 1000) {
          string refbase(chromosome_seqs[chrom] + s, 1);
          _p_svs[t].push_back(SV("INS", chrom, s, refbase, "<INS>", lc.w, 0, 0, 0, true, max(lc.l, l)));
          sa_used = true;
       } else if (l >= 2000 && l <= 50000 && lc.w >= 5) {
          string refbase(chromosome_seqs[chrom] + s, 1);
          _p_svs[t].push_back(SV("DEL", chrom, s, refbase, "<DEL>", lc.w, 0, 0, 0, true, l));
          sa_used = true;
       }
    }
    if (sa_used) continue;

    // we get the closest right clip
    int r = binary_search(rclips, 0, rclips.size() - 1, lc);
    if (r == -1) {
      continue;
    }
    auto rc = rclips[r];
    if (rc.w == 0) {
      continue;
    }

    if (abs((int)rc.p - (int)lc.p) < 1000) {
      uint s = lc.w > rc.w ? lc.p : rc.p;
      uint l = max(lc.l, rc.l);
      string refbase(chromosome_seqs[chrom] + s, 1);
      uint w = max(lc.w, rc.w);
      _p_svs[t].push_back(
          SV("INS", chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l));
    }
  }
  // Predicting deletions
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (uint i = 0; i < rclips.size(); i++) {
    const Clip &rc = rclips[i];
    int t = omp_get_thread_num();
    string chrom = rc.chrom;
    
    // For right clips, SA logic is symmetrical. We can call it here too.
    bool sa_used = false;
    if (rc.sa_has_info && rc.sa_chrom == chrom) {
       uint s = min(rc.p, rc.sa_pos);
       uint e = max(rc.p, rc.sa_pos);
       uint l = e - s;
       if (l < 1000) {
          string refbase(chromosome_seqs[chrom] + s, 1);
          _p_svs[t].push_back(SV("INS", chrom, s, refbase, "<INS>", rc.w, 0, 0, 0, true, max(rc.l, l)));
          sa_used = true;
       } else if (l >= 2000 && l <= 50000 && rc.w >= 5) {
          string refbase(chromosome_seqs[chrom] + s, 1);
          _p_svs[t].push_back(SV("DEL", chrom, s, refbase, "<DEL>", rc.w, 0, 0, 0, true, l));
          sa_used = true;
       }
    }
    if (sa_used) continue;

    // we get the closest right clip
    int l = binary_search(lclips, 0, lclips.size() - 1, rc);
    if (l == -1) {
      continue;
    }
    auto lc = lclips[l];
    if (lc.w == 0) {
      continue;
    }

    if (lc.p - rc.p >= 2000 && lc.p - rc.p <= 50000) {
      uint s = rc.p;
      uint l = lc.p - rc.p + 1;
      string refbase(chromosome_seqs[chrom] + s, 1);
      uint w = max(lc.w, rc.w);
      if (w >= 5) {
        _p_svs[t].push_back(
            SV("DEL", chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l));
      }
    }
  }
}
