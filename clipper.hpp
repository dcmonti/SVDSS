#ifndef CLIPPER_HPP
#define CLIPPER_HPP

#include <algorithm>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "interval_tree.hpp"

#include "chromosomes.hpp"
#include "sv.hpp"

using namespace std;
using namespace lib_interval_tree;

struct Clip {
  string name;
  vector<string> names;
  vector<string> sa_names;
  string chrom;
  uint p;
  uint l;
  bool starting;
  uint w;
  // Supplementary alignment info
  bool sa_has_info;
  bool primary_reverse;
  bool sa_reverse;
  string sa_chrom;
  uint sa_pos;
  uint sa_ref_len;
  uint sa_query_start;
  uint sa_query_len;
  // Weight of the winning SA group (number of reads supporting the chosen
  // breakend/event), as opposed to `w` which is the whole primary cluster.
  uint sa_w = 0;
  // Split-read support of the whole breakend, i.e. the winning group plus every
  // rival group that points at the same chromosome with the same strand pair.
  // The rivals differ from the winner only in where they place the junction, so
  // they are the same event seen with some jitter (or, for an inversion shorter
  // than SA_VOTE_POS_TOL*, its second junction). sa_w gives the geometry;
  // sa_total gives the evidence and is what the weight gate must be compared
  // against. sa_w <= sa_total <= w.
  uint sa_total = 0;
  // Inner unaligned query bases at the junction, for THIS read: the clip length
  // minus the query span its own SA segment covers. Both operands must come
  // from the same read — computing it from a cluster-level `l` (which is
  // max(l) over every merged read) against one vote-winner's sa_query_* mixes
  // provenance and yields a meaningless number.
  //
  // SIGNED on purpose, never clamped: dq > 0 means bases inserted between the
  // two segments (reported as SVINSLEN/SVINSSEQ); dq < 0 means the segments
  // overlap in query, i.e. breakpoint microhomology (reported as HOMLEN/HOMSEQ).
  // The sign is the information.
  int dq = 0;
  // Every merged read's dq, pooled, so a cluster can report the median instead
  // of one arbitrary member's value.
  vector<int> dqs;
  // The novel bases at the junction, populated only when dq > 0 (an insertion).
  // Microhomology (dq < 0) needs no storage: those bases are identical to the
  // reference at the breakpoint by definition, so HOMSEQ is read back off the
  // reference at emission time. Keeping this empty for the common
  // small-microhomology clip is what keeps the memory cost negligible.
  string ins_seq;

  Clip() { w = 0; sa_has_info = false;}

  Clip(string name_, string chrom_, const uint p_, uint l_, bool starting_,
       uint w_ = 0) {
    name = name_;
    if (!name_.empty())
      names.push_back(name_);
    chrom = chrom_;
    p = p_;
    l = l_;
    starting = starting_;
    w = w_;
    sa_has_info = false;
  }

  // Costruttore con Info SA
  Clip(string name_, string chrom_, const uint p_, uint l_, bool starting_,
       bool primary_reverse_, bool sa_reverse_, string sa_chrom_, uint sa_pos_,
       uint sa_ref_len_, uint sa_query_start_, uint sa_query_len_)
      : name(name_), chrom(chrom_), p(p_), l(l_), starting(starting_), w(1),
        sa_has_info(true), primary_reverse(primary_reverse_), sa_reverse(sa_reverse_),
        sa_chrom(sa_chrom_), sa_pos(sa_pos_), sa_ref_len(sa_ref_len_),
        sa_query_start(sa_query_start_), sa_query_len(sa_query_len_) {
    // Same-strand fallback only: the raw SA CIGAR offsets are in the SA's own
    // orientation, so this is wrong when the two segments differ in strand.
    // Clusterer::extend_alignment overwrites dq/dqs via its set_junction
    // helper, which normalises the SA span onto the primary's frame first.
    dq = (int)l_ - (int)(sa_query_start_ + sa_query_len_);
    dqs.push_back(dq);
    if (!name_.empty()) {
      names.push_back(name_);
      sa_names.push_back(name_);
    }
  }

  bool operator<(const Clip &c) const { return p < c.p; }
};

// A deletion whose pooled clip evidence is below min_cluster_weight but short
// enough (< ~15 kbp) that minimap2 may also represent it as through-reads (a
// single D op) rather than clips. Caller rescues these by counting the
// read-disjoint through-reads on the tumour BAM (Mode B).
struct ClipDelCand {
  string chrom;
  uint s;
  uint len;
  uint clip_w;          // pooled clip-side split-read support (< min_cluster_weight)
  vector<string> names; // clip read names, carried onto the emitted SV
};

class Clipper {

private:
  vector<Clip> clips;
  vector<Clip> remove_duplicates(const vector<Clip> &);
  vector<Clip> combine(const vector<Clip> &);
  vector<Clip> filter_lowcovered(const vector<Clip> &, const uint);
  vector<Clip> cluster(const vector<Clip> &, uint);
  vector<Clip> filter_tooclose_clips(
      const vector<Clip> &,
      std::unordered_map<std::string,
                         lib_interval_tree::interval_tree_t<int>> &);
  void store_clip_clusters(const vector<Clip> &lclips,
                           const vector<Clip> &rclips);

public:
  vector<vector<SV>> _p_svs;
  vector<ClipDelCand> prov_dels; // sub-threshold DELs for Mode B rescue

  Clipper(const vector<Clip> &);
  void call(int threads,
            std::unordered_map<std::string,
                               lib_interval_tree::interval_tree_t<int>> &);
};

#endif
