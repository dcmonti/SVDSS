#include "clipper.hpp"
#include "config.hpp"

namespace {
// Tolerance (bp) for grouping SA positions into the same vote during majority
// selection. Long-read split alignments at the same SV breakpoint usually
// agree to within a few tens of bp; 200 bp is a generous upper bound.
constexpr uint SA_VOTE_POS_TOL = 200;

struct SAGroup {
  string sa_chrom;
  uint sa_pos;
  uint sa_ref_len;
  uint sa_query_start;
  uint sa_query_len;
  bool primary_reverse;
  bool sa_reverse;
  uint count;
};

static void sa_vote_add(vector<SAGroup> &groups, const Clip &c) {
  if (!c.sa_has_info)
    return;
  uint w = c.w > 0 ? c.w : 1;
  for (SAGroup &g : groups) {
    if (g.sa_chrom == c.sa_chrom && g.primary_reverse == c.primary_reverse &&
        g.sa_reverse == c.sa_reverse &&
        (uint)abs((long long)g.sa_pos - (long long)c.sa_pos) <=
            SA_VOTE_POS_TOL) {
      g.count += w;
      return;
    }
  }
  groups.push_back({c.sa_chrom, c.sa_pos, c.sa_ref_len, c.sa_query_start,
                    c.sa_query_len, c.primary_reverse, c.sa_reverse, w});
}

static bool sa_vote_winner(const vector<SAGroup> &groups, SAGroup &winner) {
  if (groups.empty())
    return false;
  uint best = 0;
  for (uint i = 1; i < groups.size(); ++i) {
    if (groups[i].count > groups[best].count)
      best = i;
  }
  winner = groups[best];
  return true;
}

static void apply_sa_winner(Clip &clip, const SAGroup &w) {
  clip.sa_has_info = true;
  clip.primary_reverse = w.primary_reverse;
  clip.sa_reverse = w.sa_reverse;
  clip.sa_chrom = w.sa_chrom;
  clip.sa_pos = w.sa_pos;
  clip.sa_ref_len = w.sa_ref_len;
  clip.sa_query_start = w.sa_query_start;
  clip.sa_query_len = w.sa_query_len;
}
} // namespace

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
      vector<SAGroup> sa_votes;
      for (const Clip &c : it->second) {
        if (c.l > max_l) {
          max_l = c.l;
        }
        sa_vote_add(sa_votes, c);
      }
      Clip clip = Clip("", chrom, it->first, max_l, it->second.front().starting,
                       it->second.size());
      SAGroup winner;
      if (sa_vote_winner(sa_votes, winner))
        apply_sa_winner(clip, winner);
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
  map<uint, vector<SAGroup>> sa_votes_by_pos;
  for (const Clip &c : clips) {
    bool found = false;
    for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
         it != clusters_by_pos.end(); ++it) {
      if (it->first - r <= c.p && c.p <= it->first + r) {
        found = true;
        it->second.l = max(it->second.l, c.l);
        it->second.w += c.w;
        sa_vote_add(sa_votes_by_pos[it->first], c);
      }
    }
    if (!found) {
      Clip seed = c;
      seed.sa_has_info = false;
      clusters_by_pos[c.p] = seed;
      sa_vote_add(sa_votes_by_pos[c.p], c);
    }
  }

  for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
       it != clusters_by_pos.end(); ++it) {
    SAGroup winner;
    if (sa_vote_winner(sa_votes_by_pos[it->first], winner))
      apply_sa_winner(it->second, winner);
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
  {
    uint l_sa = 0, r_sa = 0;
    for (const Clip &c : lclips)
      if (c.sa_has_info && c.sa_query_len > 0)
        ++l_sa;
    for (const Clip &c : rclips)
      if (c.sa_has_info && c.sa_query_len > 0)
        ++r_sa;
    spdlog::info("[CLIP_SA] left: {}/{} clips with usable SA ({:.1f}%)",
                 l_sa, lclips.size(),
                 lclips.empty() ? 0.0 : 100.0 * l_sa / lclips.size());
    spdlog::info("[CLIP_SA] right: {}/{} clips with usable SA ({:.1f}%)",
                 r_sa, rclips.size(),
                 rclips.empty() ? 0.0 : 100.0 * r_sa / rclips.size());
  }
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
    if (lc.sa_has_info && lc.sa_query_len > 0) {
       uint min_sv_len = Configuration::getInstance()->min_sv_length;
       // sa_pos is 1-based (SAM spec); lc.p is 0-based (htslib) → convert
       uint sa_pos0 = lc.sa_pos > 0 ? lc.sa_pos - 1 : 0;
       if (lc.sa_chrom != chrom) {
           // BND: cross-chrom  consume clip regardless of weight
           sa_used = true;
           if (lc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + lc.p, 1);
               _p_svs[t].push_back(SV("BND", chrom, lc.p, refbase, "<BND>", lc.w, 0, 0, 0, true, 0));
           }
       } else if (lc.primary_reverse != lc.sa_reverse) {
           // INV: different strand → consume clip regardless of weight
           sa_used = true;
           uint s = min(lc.p, sa_pos0);
           uint l = max(lc.p, sa_pos0) - s;
           if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               _p_svs[t].push_back(SV("INV", chrom, s, refbase, "<INV>", lc.w, 0, 0, 0, true, l));
           }
       } else {
           // Same chrom, same strand: DEL, INS, DUP
           // lc.p is always the leftmost reference position of the primary (0-based),
           // regardless of strand. The SA for a left clip always maps to the LEFT of
           // the primary. So the reference gap is: primary_start − SA_end, for both strands.
           long long diff = (long long)lc.p - (long long)(sa_pos0 + lc.sa_ref_len);

           uint s = min(lc.p, sa_pos0);
           string refbase(chromosome_seqs[chrom] + s, 1);

           if (diff < 0) {
               // DUP
               long long jump = -diff;
               if (jump >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
                   sa_used = true;
                   _p_svs[t].push_back(SV("DUP", chrom, s, refbase, "<DUP>", lc.w, 0, 0, 0, true, jump));
               }
           } else {
               // DEL or INS
               // dR = reference gap between the two split alignments
               uint dR = (uint)diff;
               // dQ = inner unaligned query bases = clip minus SA-covered query bases
               uint sa_q = lc.sa_query_start + lc.sa_query_len;
               uint dQ = (lc.l > sa_q) ? (lc.l - sa_q) : 0;

               if (dR > dQ + min_sv_len) {
                   // DEL: reference gap exceeds query gap
                   uint l = dR - dQ;
                   if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       _p_svs[t].push_back(SV("DEL", chrom, s, refbase, "<DEL>", lc.w, 0, 0, 0, true, l));
                   }
               } else if (dQ > dR + min_sv_len) {
                   // INS: query gap exceeds reference gap
                   uint l = dQ - dR;
                   if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       _p_svs[t].push_back(SV("INS", chrom, s, refbase, "<INS>", lc.w, 0, 0, 0, true, l));
                   }
               }
               // If neither threshold met, sa_used stays false → paired-clip fallback
           }
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
      if (w >= Configuration::getInstance()->min_cluster_weight) {
        _p_svs[t].push_back(
            SV("INS", chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l));
      }
    }
  }
  // Predicting deletions
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (uint i = 0; i < rclips.size(); i++) {
    const Clip &rc = rclips[i];
    int t = omp_get_thread_num();
    string chrom = rc.chrom;
    
    // For right clips, SA logic is symmetrical.
    bool sa_used = false;
    if (rc.sa_has_info && rc.sa_query_len > 0) {
       uint min_sv_len = Configuration::getInstance()->min_sv_length;
       // sa_pos is 1-based (SAM spec); rc.p is 0-based (htslib) → convert
       uint sa_pos0 = rc.sa_pos > 0 ? rc.sa_pos - 1 : 0;
       if (rc.sa_chrom != chrom) {
           // BND: cross-chrom → consume clip regardless of weight
           sa_used = true;
           if (rc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + rc.p, 1);
               _p_svs[t].push_back(SV("BND", chrom, rc.p, refbase, "<BND>", rc.w, 0, 0, 0, true, 0));
           }
       } else if (rc.primary_reverse != rc.sa_reverse) {
           // INV: different strand → consume clip regardless of weight
           sa_used = true;
           uint target_pos = sa_pos0 + rc.sa_ref_len;
           uint s = min(rc.p, target_pos);
           uint l = max(rc.p, target_pos) - s;
           if (l >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               _p_svs[t].push_back(SV("INV", chrom, s, refbase, "<INV>", rc.w, 0, 0, 0, true, l));
           }
       } else {
           // Same chrom, same strand: DEL, INS, DUP
           // rc.p is always the rightmost reference position of the primary (0-based
           // past-end), regardless of strand. The SA for a right clip always maps to
           // the RIGHT of the primary. So the reference gap is: SA_start − primary_end,
           // for both strands.
           long long diff = (long long)sa_pos0 - (long long)rc.p;

           uint s = min(rc.p, sa_pos0);
           string refbase(chromosome_seqs[chrom] + s, 1);

           if (diff < 0) {
               // DUP
               long long jump = -diff;
               if (jump >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
                   sa_used = true;
                   _p_svs[t].push_back(SV("DUP", chrom, s, refbase, "<DUP>", rc.w, 0, 0, 0, true, jump));
               }
           } else {
               // DEL or INS
               // dR = reference gap between the two split alignments
               uint dR = (uint)diff;
               // dQ = inner unaligned query bases = clip minus SA-covered query bases
               uint sa_q = rc.sa_query_start + rc.sa_query_len;
               uint dQ = (rc.l > sa_q) ? (rc.l - sa_q) : 0;

               if (dR > dQ + min_sv_len) {
                   // DEL: reference gap exceeds query gap
                   uint l = dR - dQ;
                   if (l >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       _p_svs[t].push_back(SV("DEL", chrom, s, refbase, "<DEL>", rc.w, 0, 0, 0, true, l));
                   }
               } else if (dQ > dR + min_sv_len) {
                   // INS: query gap exceeds reference gap
                   uint l = dQ - dR;
                   if (l >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       _p_svs[t].push_back(SV("INS", chrom, s, refbase, "<INS>", rc.w, 0, 0, 0, true, l));
                   }
               }
               // If neither threshold met, sa_used stays false → paired-clip fallback
           }
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
      if (w >= Configuration::getInstance()->min_cluster_weight) {
        _p_svs[t].push_back(
            SV("DEL", chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l));
      }
    }
  }
}
