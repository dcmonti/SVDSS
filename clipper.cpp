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
  uint junction; // SA coordinate adjacent to the clip junction (vote key)
  // Per-read dq of every clip that voted for this group. The group's
  // coordinates come from the winner, so its junction length should too — and
  // as a median over the voters, not one member's value.
  vector<int> dqs;
  // Candidate inserted sequences contributed by the voters (only clips with
  // dq > 0 contribute). The representative is picked by matching length against
  // the median dq, so the emitted SVINSSEQ belongs to a read that actually
  // agrees with the reported SVINSLEN.
  vector<string> ins_seqs;
};

// Median of a per-read statistic pooled over a cluster. Median rather than
// mean or max: the distribution is tight (a real junction has every read
// agreeing to within a few bp) but carries the occasional mismapped outlier.
static int median_int(vector<int> v) {
  if (v.empty())
    return 0;
  size_t mid = v.size() / 2;
  nth_element(v.begin(), v.begin() + mid, v.end());
  return v[mid];
}

// Attach Manta/DRAGEN-style junction detail to a split-read call. These are
// reported ALONGSIDE SVLEN/END and are never netted out of them.
//
//   dq > 0 : novel bases sit between the two segments   -> SVINSLEN / SVINSSEQ
//   dq < 0 : the segments overlap in query, i.e. the shared bases align to both
//            sides                                      -> HOMLEN / HOMSEQ
//   dq = 0 : flush junction, nothing emitted (record is byte-identical to what
//            SVDSS produced before these fields existed)
//
// Homologous bases are by definition present in the reference at the junction,
// so HOMSEQ is read straight back from it instead of being carried through
// clustering — which is why Clip only ever stores inserted sequence.
static void annotate_junction(SV &sv, int dq, const string &ins_seq,
                              const char *ref_at_junction) {
  // Real breakpoint microhomology is a handful of bases; a huge value means a
  // pathological clip, not biology. There is no contig-length table in scope
  // here, so capping the copy and stopping at the terminating '\0' is what
  // keeps the read in bounds near the end of a contig.
  constexpr int MAX_HOMSEQ = 1000;
  if (dq > 0) {
    sv.ins_len = dq;
    sv.ins_seq = ins_seq;
  } else if (dq < 0) {
    sv.hom_len = -dq;
    if (ref_at_junction != nullptr) {
      int n = min(-dq, MAX_HOMSEQ);
      int k = 0;
      while (k < n && ref_at_junction[k] != '\0')
        ++k;
      sv.hom_seq.assign(ref_at_junction, (size_t)k);
    }
  }
}

// Coordinate on the SA (partner) segment that sits at the clip junction and is
// therefore consistent across reads of the same event. For a same-strand split
// the left clip's junction is the SA segment END and the right clip's is the SA
// START: a left clip's SA *start* scatters with the per-read clip length, which
// would otherwise fragment the votes of a single tandem-DUP breakend into many
// singleton groups. Opposite-strand / cross-chrom clips keep the SA start (the
// mate anchor, which does not scatter this way).
static uint sa_junction(const Clip &c) {
  if (c.sa_chrom == c.chrom && c.primary_reverse == c.sa_reverse)
    return c.starting ? (c.sa_pos + c.sa_ref_len) : c.sa_pos;
  return c.sa_pos;
}

static void sa_vote_add(vector<SAGroup> &groups, const Clip &c) {
  if (!c.sa_has_info)
    return;
  // Weight the SA vote by the number of reads that actually carry this SA (the
  // split-read support for the breakend), NOT the clip-cluster size c.w. A
  // mappability-boundary clip pile-up can stack dozens of reads that softclip at
  // the same reference position for unrelated reasons; if only one of them
  // carries the breakend's SA, using c.w would let that single split read
  // inherit the whole pile-up's weight and clear the min-cluster-weight gate,
  // emitting a high-WEIGHT phantom SV. sa_names holds exactly the SA-carrying
  // reads, so its size is the true breakend support.
  uint w = c.sa_names.size() > 0 ? (uint)c.sa_names.size() : 1;
  uint jc = sa_junction(c);
  // Pool this clip's per-read dq values into whichever group it votes for.
  const vector<int> &cdqs = c.dqs;
  for (SAGroup &g : groups) {
    if (g.sa_chrom == c.sa_chrom && g.primary_reverse == c.primary_reverse &&
        g.sa_reverse == c.sa_reverse &&
        (uint)abs((long long)g.junction - (long long)jc) <= SA_VOTE_POS_TOL) {
      g.count += w;
      g.dqs.insert(g.dqs.end(), cdqs.begin(), cdqs.end());
      if (!c.ins_seq.empty())
        g.ins_seqs.push_back(c.ins_seq);
      return;
    }
  }
  vector<string> iseqs;
  if (!c.ins_seq.empty())
    iseqs.push_back(c.ins_seq);
  groups.push_back({c.sa_chrom, c.sa_pos, c.sa_ref_len, c.sa_query_start,
                    c.sa_query_len, c.primary_reverse, c.sa_reverse, w, jc,
                    cdqs, iseqs});
}

// Pick the SA group that dictates the cluster's coordinates. `clip_chrom`/
// `clip_pos` are the cluster's own reference locus.
//
// A short event (< the 1000bp clustering radius) has both of its breakpoints
// folded into a single clip cluster: the clips physically at breakpoint A carry
// an SA pointing to B, and the clips at B point back to A. When the cluster is
// keyed at B, the A-clips' SA lands on B itself, i.e. on the cluster's own
// position, so their group predicts a zero-length event (junction == clip_pos).
// That degenerate group can out-vote the partner group (the B-clips pointing at
// A, which predict the real interval) by a single read and win, emitting
// nothing -- a 279bp inversion whose vote splits 5(self)-vs-4(partner) is the
// canonical case. A group whose junction sits within min_sv_length of the
// cluster can never yield a call, so it must not beat one that predicts a real
// interval. Fall back to the plain max-vote winner only when every group is
// degenerate (the clip is then discarded downstream anyway). Cross-chromosome
// (BND) groups are never degenerate: their junction is on another contig.
static bool sa_vote_winner(const vector<SAGroup> &groups,
                           const string &clip_chrom, uint clip_pos,
                           SAGroup &winner) {
  if (groups.empty())
    return false;
  const uint min_len = Configuration::getInstance()->min_sv_length;
  auto degenerate = [&](const SAGroup &g) {
    if (g.sa_chrom != clip_chrom)
      return false;
    long long d = (long long)g.junction - (long long)clip_pos;
    return (d < 0 ? -d : d) < (long long)min_len;
  };
  int best = -1;
  for (uint i = 0; i < groups.size(); ++i) {
    if (degenerate(groups[i]))
      continue;
    if (best < 0 || groups[i].count > groups[best].count)
      best = (int)i;
  }
  if (best < 0) {
    // every group is degenerate: keep the prior behaviour (plain max vote)
    best = 0;
    for (uint i = 1; i < groups.size(); ++i)
      if (groups[i].count > groups[best].count)
        best = (int)i;
  }
  winner = groups[best];
  return true;
}

static void apply_sa_winner(Clip &clip, const SAGroup &w,
                            const vector<SAGroup> &groups) {
  clip.sa_has_info = true;
  clip.primary_reverse = w.primary_reverse;
  clip.sa_reverse = w.sa_reverse;
  clip.sa_chrom = w.sa_chrom;
  clip.sa_pos = w.sa_pos;
  clip.sa_ref_len = w.sa_ref_len;
  clip.sa_query_start = w.sa_query_start;
  clip.sa_query_len = w.sa_query_len;
  // dq comes from the winning group too, as the median over its voters — never
  // recompute it from clip.l, which is max(l) over the whole cluster and so
  // belongs to a different read than sa_query_*.
  clip.dqs = w.dqs;
  clip.dq = median_int(w.dqs);
  // Representative inserted sequence: the one whose length matches the reported
  // median, so SVINSSEQ and SVINSLEN describe the same read. Falls back to the
  // first candidate if no voter matches exactly.
  clip.ins_seq.clear();
  if (clip.dq > 0 && !w.ins_seqs.empty()) {
    clip.ins_seq = w.ins_seqs.front();
    for (const string &s : w.ins_seqs)
      if ((int)s.size() == clip.dq) {
        clip.ins_seq = s;
        break;
      }
  }
  clip.sa_w = w.count; // reads supporting the winning event
  // Rival groups describing the same breakend must count towards the weight
  // gate, even though the winner alone dictates the coordinates. Two groups
  // describe the same breakend when they agree on the SA chromosome and on the
  // RELATIVE orientation of the two segments (primary_reverse XOR sa_reverse):
  // that XOR is what tells an inversion from a same-strand event. The two flags
  // taken separately do not, because a read sequenced from the minus strand
  // reports the very same junction with both flags mirrored — (F,T) and (T,F)
  // are the same inversion seen from the two strands, and the vote key splits
  // them into rivals for no biological reason. What is left to differ inside a
  // compatible set is where the junction sits, i.e. the jitter the vote key
  // tolerates only up to SA_VOTE_POS_TOL (and, for an inversion shorter than
  // the clip-cluster radius, its two junctions).
  const bool w_opp = (w.primary_reverse != w.sa_reverse);
  uint total = 0;
  for (const SAGroup &g : groups) {
    if (g.sa_chrom == w.sa_chrom &&
        (g.primary_reverse != g.sa_reverse) == w_opp)
      total += g.count;
  }
  clip.sa_total = total;
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

// A clip that predicts a definite reference interval for its SV type. Used to
// pool reciprocal INV/DEL evidence: an event's two breakpoints surface as two
// (or more) independent clip clusters, each below min_cluster_weight, between
// which minimap2 splits the junction reads. A read is primary on exactly one
// side, so the clusters are read-disjoint and the group weight is a plain sum.
struct IntervalCand {
  string chrom;
  uint s, e;    // predicted [start, end)
  uint w;       // split-read support to contribute (sa_w)
  bool is_left; // originating clip vector
  uint idx;     // index into that vector
};

struct IntervalGroup {
  string chrom;
  uint s, e;
  uint w;                           // summed support over the group
  vector<pair<bool, uint>> members; // (is_left, idx)
};

static inline bool within(uint a, uint b, uint tol) {
  return (uint)(a > b ? a - b : b - a) <= tol;
}

// Greedily group candidates whose predicted intervals coincide (both endpoints
// within `tol`) and sum their weights. Real events are far wider than `tol`, so
// an anchored sweep separates them cleanly. Grouping (vs a pairwise maximum)
// means an event fragmented into 3+ clusters is summed, not under-counted.
static vector<IntervalGroup> group_by_interval(const vector<IntervalCand> &cand,
                                               uint tol) {
  vector<IntervalGroup> groups;
  vector<bool> used(cand.size(), false);
  for (uint a = 0; a < cand.size(); a++) {
    if (used[a])
      continue;
    IntervalGroup g{cand[a].chrom, cand[a].s, cand[a].e, cand[a].w,
                    {{cand[a].is_left, cand[a].idx}}};
    used[a] = true;
    for (uint b = a + 1; b < cand.size(); b++) {
      if (used[b] || cand[b].chrom != g.chrom)
        continue;
      if (!within(cand[a].s, cand[b].s, tol) ||
          !within(cand[a].e, cand[b].e, tol))
        continue;
      used[b] = true;
      g.w += cand[b].w;
      g.members.push_back({cand[b].is_left, cand[b].idx});
    }
    groups.push_back(g);
  }
  return groups;
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
      if (sa_vote_winner(sa_votes, chrom, it->first, winner))
        apply_sa_winner(clip, winner, sa_votes);
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
      if (sa_vote_winner(sa_votes_by_pos[it->first], kv.first, it->first,
                         winner))
        apply_sa_winner(it->second, winner, sa_votes_by_pos[it->first]);
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

    uint dump_sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
    uint dump_dist =
        (c.p > dump_sa_pos0) ? (c.p - dump_sa_pos0) : (dump_sa_pos0 - c.p);
    uint dump_min_bnd = config->min_bnd_dist;
    bool dump_opp = (c.primary_reverse != c.sa_reverse);
    if (c.sa_chrom != c.chrom ||
        (dump_min_bnd > 0 && dump_opp && dump_dist >= dump_min_bnd)) {
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
        // Signed per-read median (see Clip::dq); clamp only for the routing
        // gate so branch selection is unchanged.
        int dq = c.dq;
        uint dQ = dq > 0 ? (uint)dq : 0u;
        dR_s = std::to_string(dR);
        dQ_s = std::to_string(dq);
        if (dR > dQ + min_sv_len) {
          branch = "DEL";
          // Mirrors the emission branch: the DEL length is the deleted
          // reference span (dR). dQ is reported in its own column.
          pred_len = std::to_string(dR);
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

  // --- Reciprocal tandem-DUP pooling ---
  // A tandem DUP surfaces as two independent clip clusters: a left clip at the
  // DUP start (its SA partner points to the DUP end) and a right clip at the DUP
  // end (its SA partner points to the DUP start). Because minimap2 splits each
  // junction read arbitrarily, the split reads are shared between the two
  // clusters, so neither side alone may reach min_cluster_weight. Detect
  // reciprocal DUP-voting clusters and pool their SA weights so the combined
  // evidence can clear the gate. Both loops then emit the same DUP (identical
  // coordinates), collapsed by the RO>=0.9 dedup in Caller::run().
  vector<uint> lclip_pool(lclips.size(), 0);
  vector<uint> rclip_pool(rclips.size(), 0);
  {
    const uint RECIP_TOL = 2000; // breakpoint match tolerance (bp)
    // A same-chrom/same-strand clip whose geometry yields a DUP (diff<0). Sets
    // `nearp` (the clip-side breakend) and `partner` (the SA-side breakend).
    auto dup_anchors = [](const Clip &c, bool is_left, uint &nearp,
                          uint &partner) -> bool {
      if (!c.sa_has_info || c.sa_query_len == 0) return false;
      if (c.sa_chrom != c.chrom) return false;
      if (c.primary_reverse != c.sa_reverse) return false;
      uint sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
      long long diff = is_left
          ? (long long)c.p - (long long)(sa_pos0 + c.sa_ref_len)
          : (long long)sa_pos0 - (long long)c.p;
      if (diff >= 0) return false;
      nearp = c.p;
      partner = is_left ? (sa_pos0 + c.sa_ref_len) : sa_pos0;
      return true;
    };
    for (uint i = 0; i < lclips.size(); i++) {
      uint lnear, lpart;
      if (!dup_anchors(lclips[i], true, lnear, lpart)) continue;
      for (uint j = 0; j < rclips.size(); j++) {
        if (rclips[j].chrom != lclips[i].chrom) continue;
        uint rnear, rpart;
        if (!dup_anchors(rclips[j], false, rnear, rpart)) continue;
        // reciprocal: left clip-side == right partner (DUP start) AND
        //             left partner == right clip-side (DUP end).
        if ((uint)abs((long long)lnear - (long long)rpart) <= RECIP_TOL &&
            (uint)abs((long long)lpart - (long long)rnear) <= RECIP_TOL) {
          uint pooled = lclips[i].sa_w + rclips[j].sa_w;
          if (pooled > lclip_pool[i]) lclip_pool[i] = pooled;
          if (pooled > rclip_pool[j]) rclip_pool[j] = pooled;
        }
      }
    }
  }

  // Reciprocal INVERSION + DELETION pooling. Both events have two breakpoints
  // between which minimap2 splits the junction reads, so neither clip cluster
  // may reach min_cluster_weight alone — the failure the tandem-DUP pooling
  // above fixes for duplications. Each cluster independently predicts the SAME
  // reference interval (via the emission formulas below), so we pair by
  // predicted interval and pool the split-read weights. DUP is handled
  // separately above: its clip coordinates are imprecise, so it matches
  // breakends reciprocally instead of predicting a clean interval.
  vector<uint> linv_pool(lclips.size(), 0), rinv_pool(rclips.size(), 0);
  vector<uint> ldel_pool(lclips.size(), 0), rdel_pool(rclips.size(), 0);
  {
    const uint min_sv_len = Configuration::getInstance()->min_sv_length;

    // Predicted inversion interval (opposite strand): left clip -> sa_pos,
    // right clip -> sa_pos + sa_ref_len (mirrors the emission code).
    auto inv_pred = [](const Clip &c, bool is_left, uint &s, uint &e) -> bool {
      if (!c.sa_has_info || c.sa_query_len == 0) return false;
      if (c.sa_chrom != c.chrom) return false;
      if (c.primary_reverse == c.sa_reverse) return false; // opposite strand
      uint sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
      uint target = is_left ? sa_pos0 : (sa_pos0 + c.sa_ref_len);
      s = min(c.p, target);
      e = max(c.p, target);
      return e > s; // a degenerate (zero-length) prediction carries no interval
    };

    // Predicted deletion interval (same strand, diff > 0). The start is the
    // breakpoint each emission branch uses: sa_pos0 + sa_ref_len for a left
    // clip, the primary end (== c.p, since diff > 0) for a right clip.
    auto del_pred = [&](const Clip &c, bool is_left, uint &s, uint &e) -> bool {
      if (!c.sa_has_info || c.sa_query_len == 0) return false;
      if (c.sa_chrom != c.chrom) return false;
      if (c.primary_reverse != c.sa_reverse) return false; // same strand
      uint sa_pos0 = c.sa_pos > 0 ? c.sa_pos - 1 : 0;
      long long diff =
          is_left ? (long long)c.p - (long long)(sa_pos0 + c.sa_ref_len)
                  : (long long)sa_pos0 - (long long)c.p;
      if (diff <= 0) return false; // DUP side (diff < 0) handled elsewhere
      uint dR = (uint)diff;
      uint dQ = c.dq > 0 ? (uint)c.dq : 0u; // clamped: routing gate only
      if (dR <= dQ + min_sv_len) return false; // INS side / too short for a DEL
      uint start = is_left ? (sa_pos0 + c.sa_ref_len) : c.p;
      s = start;
      // Must mirror the emission branches exactly: the interval spans the
      // deleted reference (dR), not the net dR - dQ. A short end here shrinks
      // the pooling window and can stop the two reciprocal breakpoints of one
      // deletion from pooling at all.
      e = start + dR;
      return e > s;
    };

    // Flatten both clip sides into one candidate list so L/L, R/R and L/R pairs
    // are all covered.
    auto build = [&](auto pred) {
      vector<IntervalCand> cand;
      for (uint i = 0; i < lclips.size(); i++) {
        uint s, e;
        if (pred(lclips[i], true, s, e))
          cand.push_back({lclips[i].chrom, s, e, lclips[i].sa_w, true, i});
      }
      for (uint j = 0; j < rclips.size(); j++) {
        uint s, e;
        if (pred(rclips[j], false, s, e))
          cand.push_back({rclips[j].chrom, s, e, rclips[j].sa_w, false, j});
      }
      return cand;
    };
    auto apply_pool = [&](const vector<IntervalGroup> &gs, vector<uint> &lp,
                          vector<uint> &rp) {
      for (const IntervalGroup &g : gs)
        for (const auto &m : g.members) {
          uint &slot = m.first ? lp[m.second] : rp[m.second];
          if (g.w > slot) slot = g.w;
        }
    };

    apply_pool(group_by_interval(build(inv_pred), 100), linv_pool, rinv_pool);

    vector<IntervalGroup> del_groups = group_by_interval(build(del_pred), 500);
    apply_pool(del_groups, ldel_pool, rdel_pool);

    // Mode B candidates: a DEL group that stays below the gate on clips alone
    // AND is short enough that minimap2 may also have emitted through-reads
    // (< 15 kbp; it does not open D ops much beyond ~10 kbp). Caller adds the
    // read-disjoint through-read count on the tumour BAM.
    for (const IntervalGroup &g : del_groups) {
      uint len = g.e - g.s;
      if (g.w >= min_cw) continue;                    // already clears on clips
      if (len < min_sv_len || len >= 15000) continue; // Mode B window
      ClipDelCand pc;
      pc.chrom = g.chrom;
      pc.s = g.s;
      pc.len = len;
      pc.clip_w = g.w;
      for (const auto &m : g.members) {
        const Clip &c = m.first ? lclips[m.second] : rclips[m.second];
        append_names(pc.names, c.names);
      }
      prov_dels.push_back(std::move(pc));
    }
  }

  // Predicting insertions
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (uint i = 0; i < lclips.size(); i++) {
    const Clip &lc = lclips[i];
    int t = omp_get_thread_num();
    string chrom = lc.chrom;

    bool sa_used = false;
    if (lc.sa_has_info && lc.sa_query_len > 0) {
       // Weight gate on the split-read support of the breakend: if too few reads
       // carry a supplementary alignment across this junction, discard the clip
       // (and do not fall through to the paired-clip fallback, since this is an
       // SA-carrying clip).
       // The support is sa_total, not the winning group's sa_w: rival groups on
       // the same chromosome and strand pair are the same breakend placed a few
       // bp apart (the vote only merges junctions within SA_VOTE_POS_TOL), so
       // gating on the winner alone throws away real split reads and drops
       // events that have plenty of them — a 279 bp inversion, whose two
       // junctions land in two rival groups, is a clean example.
       // eff_w folds in reciprocal pooling: for a DUP-voting clip that found a
       // reciprocal partner, lclip_pool[i] holds the summed (left+right)
       // split-read weight; linv_pool[i] does the same for an inversion whose
       // two junction clusters predict the same interval. Both are 0 when no
       // partner was found, so eff_w falls back to sa_total.
       uint eff_w = max(lc.sa_total,
                        max(lclip_pool[i], max(linv_pool[i], ldel_pool[i])));
       if (eff_w < Configuration::getInstance()->min_cluster_weight)
         continue;
       uint min_sv_len = Configuration::getInstance()->min_sv_length;
       // sa_pos is 1-based (SAM spec); lc.p is 0-based (htslib) → convert
       uint sa_pos0 = lc.sa_pos > 0 ? lc.sa_pos - 1 : 0;
       // Same-chrom split alignments farther apart than min_bnd_dist are a
       // long-range (intra-chromosomal) translocation: emit a BND breakend like
       // cross-chrom, not one giant INV/DEL/DUP spanning the whole gap.
       uint intra_dist = (lc.p > sa_pos0) ? (lc.p - sa_pos0) : (sa_pos0 - lc.p);
       uint min_bnd_dist = Configuration::getInstance()->min_bnd_dist;
       // Only OPPOSITE-strand long-range junctions become BND: a huge same-strand
       // event is a genuine DEL/DUP (a 10Mbp deletion must stay a DEL), while a
       // huge opposite-strand "inversion" is really a translocation junction.
       bool opp_strand = (lc.primary_reverse != lc.sa_reverse);
       bool as_bnd = (lc.sa_chrom != chrom) ||
                     (min_bnd_dist > 0 && opp_strand && intra_dist >= min_bnd_dist);
       if (as_bnd) {
           // BND: cross-chrom or long-range intra-chrom translocation.
           // Consume clip regardless of weight.
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
           // Gate on eff_w, not lc.w: when the inversion's two junction clusters
           // each hold part of the split-read support, only the pooled weight
           // reflects the true evidence. Both loops then emit the same INV
           // (identical coordinates), collapsed by the RO>=0.9 dedup in
           // Caller::run(), which also picks the higher-WEIGHT copy.
           uint inv_w = max(lc.w, linv_pool[i]);
           if (l >= min_sv_len && inv_w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               SV sv = SV("INV", chrom, s, refbase, "<INV>", inv_w, 0, 0, 0, true, l);
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

           if (diff < 0) {
               // DUP — duplicated interval spans the primary/SA overlap. Geometry
               // differs from DEL/INS; left unchanged (its own coordinate handling).
               uint s = min(lc.p, sa_pos0);
               string refbase(chromosome_seqs[chrom] + s, 1);
               long long jump = -diff;
               if (jump >= min_sv_len && eff_w >= Configuration::getInstance()->min_cluster_weight) {
                   sa_used = true;
                   SV sv = SV("DUP", chrom, s, refbase, "<DUP>", eff_w, 0, 0, 0, true, jump);
                   sv.add_reads(lc.names);
                   sv.add_sa_reads(lc.sa_names);
                   _p_svs[t].push_back(sv);
               }
           } else {
               // DEL or INS. The breakpoint adjacent to the junction is the SA
               // segment END (sa_pos0 + sa_ref_len), NOT the SA start: the SA start
               // scatters with the per-read clip length, so using it places the
               // event ~sa_ref_len bp upstream of where the reads actually clip,
               // producing a shifted phantom duplicate of the right-clip call.
               uint j = sa_pos0 + lc.sa_ref_len;
               string refbase(chromosome_seqs[chrom] + j, 1);
               // dR = reference gap between the two split alignments
               uint dR = (uint)diff;
               // dq = inner unaligned query bases at the junction, as the median
               // over the reads that voted for this SA group (see Clip::dq).
               // Signed: >0 inserted bases, <0 breakpoint microhomology.
               int dq = lc.dq;
               // Clamped copy for the routing gate: a microhomology junction
               // (dq < 0) must not make the gate MORE permissive than a flush
               // one, so negatives are treated as 0 exactly as before.
               //
               // Note this does NOT keep branch selection identical to the old
               // code — it cannot. dq's *value* is now correct rather than
               // mixed-provenance garbage, so events do move between the DEL,
               // INS and NONE branches. That is the point: chr14:37,300,397
               // was routed to INS by a bogus dQ=3038 and is a DEL(1625) once
               // dq is computed honestly (1271).
               uint dQ = dq > 0 ? (uint)dq : 0u;

               if (dR > dQ + min_sv_len) {
                   // DEL: reference gap exceeds query gap.
                   // The length handed to SV() must be the DELETED REFERENCE SPAN
                   // (dR), not the net allele-length change (dR - dQ). SV() derives
                   // e = s + l, and j + dR == lc.p exactly — the observed distal
                   // breakpoint — so netting dQ out here places END dQ bp short of
                   // where the reads actually resume. Per VCF 4.4, SVLEN for a
                   // symbolic <DEL> is the number of deleted reference bases; any
                   // bases inserted at the junction are reported separately via
                   // SVINSLEN/SVINSSEQ (see dQ handling below).
                   uint l = dR;
                   // Gate on the pooled weight: when the deletion's two junction
                   // clusters each hold part of the split-read support, only the
                   // pooled weight reflects the true evidence (see INV/DUP).
                   uint del_w = max(lc.w, ldel_pool[i]);
                   if (l >= min_sv_len && del_w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       SV sv = SV("DEL", chrom, j, refbase, "<DEL>", del_w, 0, 0, 0, true, l);
                       annotate_junction(sv, dq, lc.ins_seq,
                                         chromosome_seqs[chrom] + j);
                       sv.add_reads(lc.names);
                       sv.add_sa_reads(lc.sa_names);
                       _p_svs[t].push_back(sv);
                   }
               } else if (dQ > dR + min_sv_len) {
                   // INS: query gap exceeds reference gap
                   uint l = dQ - dR;
                   if (l >= min_sv_len && lc.w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       SV sv = SV("INS", chrom, j, refbase, "<INS>", lc.w, 0, 0, 0, true, l);
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
       // Weight gate on the breakend's split-read support (see left-clip loop):
       // rc.sa_total, i.e. the winning SA group plus its same-chrom same-strand
       // rivals, folded together with reciprocal DUP and INV pooling.
       uint eff_w = max(rc.sa_total,
                        max(rclip_pool[i], max(rinv_pool[i], rdel_pool[i])));
       if (eff_w < Configuration::getInstance()->min_cluster_weight)
         continue;
       uint min_sv_len = Configuration::getInstance()->min_sv_length;
       // sa_pos is 1-based (SAM spec); rc.p is 0-based (htslib) → convert
       uint sa_pos0 = rc.sa_pos > 0 ? rc.sa_pos - 1 : 0;
       // Same-chrom split alignments farther apart than min_bnd_dist are a
       // long-range (intra-chromosomal) translocation: emit a BND breakend like
       // cross-chrom, not one giant INV/DEL/DUP spanning the whole gap.
       uint intra_dist = (rc.p > sa_pos0) ? (rc.p - sa_pos0) : (sa_pos0 - rc.p);
       uint min_bnd_dist = Configuration::getInstance()->min_bnd_dist;
       // Only OPPOSITE-strand long-range junctions become BND: a huge same-strand
       // event is a genuine DEL/DUP (a 10Mbp deletion must stay a DEL), while a
       // huge opposite-strand "inversion" is really a translocation junction.
       bool opp_strand = (rc.primary_reverse != rc.sa_reverse);
       bool as_bnd = (rc.sa_chrom != chrom) ||
                     (min_bnd_dist > 0 && opp_strand && intra_dist >= min_bnd_dist);
       if (as_bnd) {
           // BND: cross-chrom or long-range intra-chrom translocation.
           // Consume clip regardless of weight.
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
           // Gate on the pooled weight (see left-clip loop).
           uint inv_w = max(rc.w, rinv_pool[i]);
           if (l >= min_sv_len && inv_w >= Configuration::getInstance()->min_cluster_weight) {
               string refbase(chromosome_seqs[chrom] + s, 1);
               SV sv = SV("INV", chrom, s, refbase, "<INV>", inv_w, 0, 0, 0, true, l);
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
               if (jump >= min_sv_len && eff_w >= Configuration::getInstance()->min_cluster_weight) {
                   sa_used = true;
                   SV sv = SV("DUP", chrom, s, refbase, "<DUP>", eff_w, 0, 0, 0, true, jump);
                   sv.add_reads(rc.names);
                   sv.add_sa_reads(rc.sa_names);
                   _p_svs[t].push_back(sv);
               }
           } else {
               // DEL or INS
               // dR = reference gap between the two split alignments
               uint dR = (uint)diff;
               // dq = per-read median inner unaligned query bases (see Clip::dq).
               int dq = rc.dq;
               // Clamped for the routing gate only; see the left-clip loop for
               // why this does not preserve the old branch selection.
               uint dQ = dq > 0 ? (uint)dq : 0u;

               if (dR > dQ + min_sv_len) {
                   // DEL: reference gap exceeds query gap. l is the DELETED
                   // REFERENCE SPAN (dR), not the net dR - dQ — see the left-clip
                   // loop for the full rationale. Here s + dR == sa_pos0 exactly,
                   // the observed distal breakpoint.
                   uint l = dR;
                   // Gate on the pooled weight (see left-clip loop).
                   uint del_w = max(rc.w, rdel_pool[i]);
                   if (l >= min_sv_len && del_w >= Configuration::getInstance()->min_cluster_weight) {
                       sa_used = true;
                       SV sv = SV("DEL", chrom, s, refbase, "<DEL>", del_w, 0, 0, 0, true, l);
                       annotate_junction(sv, dq, rc.ins_seq,
                                         chromosome_seqs[chrom] + s);
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
