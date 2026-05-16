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
    if (!name_.empty()) {
      names.push_back(name_);
      sa_names.push_back(name_);
    }
  }

  bool operator<(const Clip &c) const { return p < c.p; }
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

  Clipper(const vector<Clip> &);
  void call(int threads,
            std::unordered_map<std::string,
                               lib_interval_tree::interval_tree_t<int>> &);
};

#endif
