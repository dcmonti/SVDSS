#include "clipper.hpp"
#include "config.hpp"

#include <fstream>

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

static void append_names(vector<string> &out, const vector<string> &in) {
  for (const string &name : in) {
    if (!name.empty())
      out.push_back(name);
  }
}

static void append_name(vector<string> &out, const string &name) {
  if (!name.empty())
    out.push_back(name);
}

static string join_names(const vector<string> &v) {
  if (v.empty())
    return ".";
  string out;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i)
      out.push_back(',');
    out += v[i];
  }
  return out;
}

static string strands_to_string(bool primary_reverse, bool sa_reverse) {
  string s = "?/?";
  s[0] = primary_reverse ? '-' : '+';
  s[2] = sa_reverse ? '-' : '+';
  return s;
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
      vector<string> names;
      vector<string> sa_names;
      for (const Clip &c : it->second) {
        if (c.l > max_l) {
          max_l = c.l;
        }
        if (!c.names.empty())
          append_names(names, c.names);
        else
          append_name(names, c.name);
        if (c.sa_has_info) {
          if (!c.sa_names.empty())
            append_names(sa_names, c.sa_names);
          else
            append_name(sa_names, c.name);
        }
        sa_vote_add(sa_votes, c);
      }
      Clip clip = Clip("", chrom, it->first, max_l, it->second.front().starting,
                       it->second.size());
      clip.names = names;
      clip.sa_names = sa_names;
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

// Cluster clips by proximity. Clustering is performed independently for each
// chromosome so that nearby positions on different chroms cannot be merged.
// TODO: this might be too slow
vector<Clip> Clipper::cluster(const vector<Clip> &clips, uint r) {
  unordered_map<string, vector<Clip>> by_chrom;
  for (const Clip &c : clips)
    by_chrom[c.chrom].push_back(c);

  vector<Clip> clusters;
  for (auto &kv : by_chrom) {
    map<uint, Clip> clusters_by_pos;
    map<uint, vector<SAGroup>> sa_votes_by_pos;
    for (const Clip &c : kv.second) {
      bool found = false;
      for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
           it != clusters_by_pos.end(); ++it) {
        if (it->first - r <= c.p && c.p <= it->first + r) {
          found = true;
          it->second.l = max(it->second.l, c.l);
          it->second.w += c.w;
          if (!c.names.empty())
            append_names(it->second.names, c.names);
          else
            append_name(it->second.names, c.name);
          if (c.sa_has_info) {
            if (!c.sa_names.empty())
              append_names(it->second.sa_names, c.sa_names);
            else
              append_name(it->second.sa_names, c.name);
          }
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
  }
  return clusters;
}

// Dump the final, fully-filtered/clustered left and right clip sets to a TSV
// file. For SA-carrying clips also pre-compute the branch (BND/INV/DUP/DEL/INS
// or NONE) and predicted length using the same formulas as Path 1 in
// Clipper::call, so the file is enough to debug routing decisions without
// rerunning with debug logs. No weight thresholds applied here: all final
// clusters are written regardless of min_cluster_weight.
void Clipper::store_clip_clusters(const vector<Clip> &lclips,
                                  const vector<Clip> &rclips) {
  Configuration *config = Configuration::getInstance();
  if (config->clip_clusters.empty())
    return;
  ofstream f;
  f.open(config->clip_clusters);
  f << "#SVDSS clip clusters\n";
  f << "#fields=chrom:p_1based\tside\tw\tl_max\tn_reads\tn_sa_reads"
       "\tstrands\tsa_chrom:sa_pos_1based\tsa_ref_len\tsa_query_start"
       "\tsa_query_len\tdiff\tdR\tdQ\tbranch\tpred_len\treads\tsa_reads\n";

  uint min_sv_len = config->min_sv_length;

  auto write_one = [&](const Clip &c, bool is_left) {
    f << c.chrom << ":" << (c.p + 1) << "\t" << (is_left ? "L" : "R") << "\t"
      << c.w << "\t" << c.l << "\t" << c.names.size() << "\t"
      << c.sa_names.size();

    if (!c.sa_has_info || c.sa_query_len == 0) {
      // No usable SA info: paired-clip fallback territory.
      for (int k = 0; k < 9; ++k)
        f << "\t.";
      f << "\t" << join_names(c.names) << "\t" << join_names(c.sa_names)
        << "\n";
      return;
    }

    f << "\t" << strands_to_string(c.primary_reverse, c.sa_reverse) << "\t"
      << c.sa_chrom << ":" << c.sa_pos << "\t" << c.sa_ref_len << "\t"
      << c.sa_query_start << "\t" << c.sa_query_len;

    string branch;
    string pred_len = ".";
    string diff_s = ".", dR_s = ".", dQ_s = ".";

    if (c.sa_chrom != c.chrom) {
      branch = "BND";
    } else if (c.primary_reverse != c.sa_reverse) {
      branch = "INV";
      uint sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
      uint inv_len = 0;
      if (is_left) {
        uint s = std::min(c.p, sa_pos0);
        inv_len = std::max(c.p, sa_pos0) - s;
      } else {
        uint target_pos = sa_pos0 + c.sa_ref_len;
        uint s = std::min(c.p, target_pos);
        inv_len = std::max(c.p, target_pos) - s;
      }
      pred_len = std::to_string(inv_len);
    } else {
      uint sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
      long long diff =
          is_left ? (long long)c.p - (long long)(sa_pos0 + c.sa_ref_len)
                  : (long long)sa_pos0 - (long long)c.p;
      diff_s = std::to_string(diff);
      if (diff < 0) {
        branch = "DUP";
        pred_len = std::to_string(-diff);
      } else {
        uint dR = (uint)diff;
        uint sa_q = c.sa_query_start + c.sa_query_len;
        uint dQ = (c.l > sa_q) ? (c.l - sa_q) : 0;
        dR_s = std::to_string(dR);
        dQ_s = std::to_string(dQ);
        if (dR > dQ + min_sv_len) {
          branch = "DEL";
          pred_len = std::to_string(dR - dQ);
        } else if (dQ > dR + min_sv_len) {
          branch = "INS";
          pred_len = std::to_string(dQ - dR);
        } else {
          // Same chrom/strand but neither threshold met → falls through to
          // paired-clip fallback if enabled.
          branch = "NONE";
        }
      }
    }

    f << "\t" << diff_s << "\t" << dR_s << "\t" << dQ_s << "\t" << branch
      << "\t" << pred_len << "\t" << join_names(c.names) << "\t"
      << join_names(c.sa_names) << "\n";
  };

  for (const Clip &c : lclips)
    write_one(c, true);
  for (const Clip &c : rclips)
    write_one(c, false);
  f.close();
}

vector<Clip> Clipper::filter_tooclose_clips(
    const vector<Clip> &clips,
    unordered_map<string, interval_tree_t<int>> &vartrees) {
  vector<Clip> fclips;
  for (const Clip &c : clips) {
    auto it = vartrees.find(c.chrom);
    if (it == vartrees.end() ||
        it->second.overlap_find({(int)c.p, (int)c.p + 1}) ==
            end(it->second)) {
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

void Clipper::call(int threads,
                   unordered_map<string, interval_tree_t<int>> &vartrees) {
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
      rclips = filter_tooclose_clips(rclips, vartrees);
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
      lclips = filter_tooclose_clips(lclips, vartrees);
      spdlog::info("[CLIP_FILTER][LEFT] after filter_tooclose_clips: {}", lclips.size());
      lclips = cluster(lclips, 1000);
      spdlog::info("[CLIP_FILTER][LEFT] after cluster(1000): {}", lclips.size());
      sort(lclips.begin(), lclips.end());
    }
  }
  // Per-chromosome views, used by Path 2 paired-clip lookup so that
  // binary_search on `p` can never cross chromosome boundaries. Each per-chrom
  // vector inherits the global sort order.
  unordered_map<string, vector<Clip>> lclips_by_chrom;
  unordered_map<string, vector<Clip>> rclips_by_chrom;
  for (const Clip &c : lclips)
    lclips_by_chrom[c.chrom].push_back(c);
  for (const Clip &c : rclips)
    rclips_by_chrom[c.chrom].push_back(c);
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
  if (!Configuration::getInstance()->clip_clusters.empty()) {
    spdlog::info("Storing clip clusters to {}",
                 Configuration::getInstance()->clip_clusters);
    store_clip_clusters(lclips, rclips);
  }
  _p_svs.resize(threads);
  // Aggregate clip-cluster statistics (info-level summary at the end).
  uint min_cw = Configuration::getInstance()->min_cluster_weight;
  size_t l_above_thr = 0, r_above_thr = 0; // clip clusters with w >= threshold
  size_t l_with_call = 0, r_with_call = 0; // clip clusters producing >=1 SV
  for (const Clip &c : lclips)
    if (c.w >= min_cw)
      ++l_above_thr;
  for (const Clip &c : rclips)
    if (c.w >= min_cw)
      ++r_above_thr;
  if (lclips.empty() || rclips.empty()) {
    spdlog::info(
        "[CLIP_STATS] left clusters: total={} above_weight(>={})={} produced_call=0 | "
        "right clusters: total={} above_weight(>={})={} produced_call=0",
        lclips.size(), min_cw, l_above_thr, rclips.size(), min_cw, r_above_thr);
    return;
  }
  // Each clip produces at most one SV in these loops, so the number of clip
  // clusters that produced a call equals the total SVs added by the loop.
  auto total_svs = [&]() {
    size_t n = 0;
    for (int t = 0; t < threads; ++t)
      n += _p_svs[t].size();
    return n;
  };
  size_t svs_before_left = total_svs();
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
               // Mate junction position on the partner chromosome. For a left
               // clip the SA end adjacent to the junction is sa_pos+sa_ref_len
               // when both map on the same strand, otherwise sa_pos (1-based).
               bool same_strand = (lc.primary_reverse == lc.sa_reverse);
               uint mate = same_strand ? (lc.sa_pos + lc.sa_ref_len) : lc.sa_pos;
               // Left clip: primary lies to the RIGHT of the breakend, so the
               // mate piece precedes refbase. ']' keeps the mate forward (same
               // strand), '[' takes it reverse-complemented (opposite strand).
               string br = same_strand ? "]" : "[";
               string alt = br + lc.sa_chrom + ":" + to_string(mate) + br + refbase;
               SV sv = SV("BND", chrom, lc.p, refbase, alt, lc.w, 0, 0, 0, true, 0);
               sv.add_reads(lc.names);
               sv.add_sa_reads(lc.sa_names);
               _p_svs[t].push_back(sv);
           }
       } else if (lc.primary_reverse != lc.sa_reverse) {
           // INV: different strand → consume clip regardless of weight
           sa_used = true;
           uint s = min(lc.p, sa_pos0);
           uint l = max(lc.p, sa_pos0) - s;
           if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               SV sv = SV("INV", chrom, s, refbase, "<INV>", lc.w, 0, 0, 0, true, l);
               sv.add_reads(lc.names);
               sv.add_sa_reads(lc.sa_names);
               _p_svs[t].push_back(sv);
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
                   SV sv = SV("DUP", chrom, s, refbase, "<DUP>", lc.w, 0, 0, 0, true, jump);
                   sv.add_reads(lc.names);
                   sv.add_sa_reads(lc.sa_names);
                   _p_svs[t].push_back(sv);
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
                       SV sv = SV("DEL", chrom, s, refbase, "<DEL>", lc.w, 0, 0, 0, true, l);
                       sv.add_reads(lc.names);
                       sv.add_sa_reads(lc.sa_names);
                       _p_svs[t].push_back(sv);
                   }
               } else if (dQ > dR + min_sv_len) {
                   // INS: query gap exceeds reference gap
                   uint l = dQ - dR;
                   if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       SV sv = SV("INS", chrom, s, refbase, "<INS>", lc.w, 0, 0, 0, true, l);
                       sv.add_reads(lc.names);
                       sv.add_sa_reads(lc.sa_names);
                       _p_svs[t].push_back(sv);
                   }
               }
               // If neither threshold met, sa_used stays false → paired-clip fallback
           }
       }
    }
    if (sa_used) continue;

    // If SA evidence exists for this clip, skip paired-clip fallback to avoid
    // emitting conflicting calls on the same cluster.
    if (lc.sa_has_info && lc.sa_query_len > 0)
      continue;

    if (!Configuration::getInstance()->clipped_fallback)
      continue;

    // we get the closest right clip on the same chromosome
    auto rit = rclips_by_chrom.find(lc.chrom);
    if (rit == rclips_by_chrom.end() || rit->second.empty())
      continue;
    int r = binary_search(rit->second, 0, rit->second.size() - 1, lc);
    if (r == -1) {
      continue;
    }
    auto rc = rit->second[r];
    if (rc.w == 0) {
      continue;
    }

    if (abs((int)rc.p - (int)lc.p) < 1000) {
      uint s = lc.w > rc.w ? lc.p : rc.p;
      uint l = max(lc.l, rc.l);
      string refbase(chromosome_seqs[chrom] + s, 1);
      uint w = max(lc.w, rc.w);
      if (w >= Configuration::getInstance()->min_cluster_weight) {
        vector<string> pair_names;
        append_names(pair_names, lc.names);
        append_names(pair_names, rc.names);
        SV sv = SV("INS", chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l);
        sv.add_reads(pair_names);
        _p_svs[t].push_back(sv);
      }
    }
  }
  l_with_call = total_svs() - svs_before_left;
  size_t svs_before_right = total_svs();
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
               // Mate junction position: symmetrical to the left-clip case. For
               // a right clip the adjacent SA end is sa_pos when both map on the
               // same strand, otherwise sa_pos+sa_ref_len (1-based).
               bool same_strand = (rc.primary_reverse == rc.sa_reverse);
               uint mate = same_strand ? rc.sa_pos : (rc.sa_pos + rc.sa_ref_len);
               // Right clip: primary lies to the LEFT of the breakend, so the
               // mate piece follows refbase. '[' keeps the mate forward (same
               // strand), ']' takes it reverse-complemented (opposite strand).
               string br = same_strand ? "[" : "]";
               string alt = refbase + br + rc.sa_chrom + ":" + to_string(mate) + br;
               SV sv = SV("BND", chrom, rc.p, refbase, alt, rc.w, 0, 0, 0, true, 0);
               sv.add_reads(rc.names);
               sv.add_sa_reads(rc.sa_names);
               _p_svs[t].push_back(sv);
           }
       } else if (rc.primary_reverse != rc.sa_reverse) {
           // INV: different strand → consume clip regardless of weight
           sa_used = true;
           uint target_pos = sa_pos0 + rc.sa_ref_len;
           uint s = min(rc.p, target_pos);
           uint l = max(rc.p, target_pos) - s;
           if (l >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               SV sv = SV("INV", chrom, s, refbase, "<INV>", rc.w, 0, 0, 0, true, l);
               sv.add_reads(rc.names);
               sv.add_sa_reads(rc.sa_names);
               _p_svs[t].push_back(sv);
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
                   SV sv = SV("DUP", chrom, s, refbase, "<DUP>", rc.w, 0, 0, 0, true, jump);
                   sv.add_reads(rc.names);
                   sv.add_sa_reads(rc.sa_names);
                   _p_svs[t].push_back(sv);
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
                       SV sv = SV("DEL", chrom, s, refbase, "<DEL>", rc.w, 0, 0, 0, true, l);
                       sv.add_reads(rc.names);
                       sv.add_sa_reads(rc.sa_names);
                       _p_svs[t].push_back(sv);
                   }
               } else if (dQ > dR + min_sv_len) {
                   // INS: query gap exceeds reference gap
                   uint l = dQ - dR;
                   if (l >= min_sv_len && rc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       SV sv = SV("INS", chrom, s, refbase, "<INS>", rc.w, 0, 0, 0, true, l);
                       sv.add_reads(rc.names);
                       sv.add_sa_reads(rc.sa_names);
                       _p_svs[t].push_back(sv);
                   }
               }
               // If neither threshold met, sa_used stays false → paired-clip fallback
           }
       }
    }
    if (sa_used) continue;

    // If SA evidence exists for this clip, skip paired-clip fallback to avoid
    // emitting conflicting calls on the same cluster.
    if (rc.sa_has_info && rc.sa_query_len > 0)
      continue;

    if (!Configuration::getInstance()->clipped_fallback)
      continue;

    // we get the closest left clip on the same chromosome
    auto lit = lclips_by_chrom.find(rc.chrom);
    if (lit == lclips_by_chrom.end() || lit->second.empty())
      continue;
    int l = binary_search(lit->second, 0, lit->second.size() - 1, rc);
    if (l == -1) {
      continue;
    }
    auto lc = lit->second[l];
    if (lc.w == 0) {
      continue;
    }

    if (lc.p - rc.p >= 2000 && lc.p - rc.p <= 50000) {
      uint s = rc.p;
      uint l = lc.p - rc.p + 1;
      string refbase(chromosome_seqs[chrom] + s, 1);
      uint w = max(lc.w, rc.w);
      if (w >= Configuration::getInstance()->min_cluster_weight) {
        vector<string> pair_names;
        append_names(pair_names, lc.names);
        append_names(pair_names, rc.names);
        SV sv = SV("DEL", chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l);
        sv.add_reads(pair_names);
        _p_svs[t].push_back(sv);
      }
    }
  }
  r_with_call = total_svs() - svs_before_right;
  spdlog::info(
      "[CLIP_STATS] left clusters: total={} above_weight(>={})={} produced_call={} | "
      "right clusters: total={} above_weight(>={})={} produced_call={}",
      lclips.size(), min_cw, l_above_thr, l_with_call, rclips.size(), min_cw,
      r_above_thr, r_with_call);
}
