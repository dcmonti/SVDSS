#ifndef SV_HPP
#define SV_HPP

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class SV {
public:
  string type;
  string chrom;
  string idx;
  int s;
  int e;
  string refall;
  string altall;
  uint w;
  int cov;
  int cov0;
  int cov1;
  int cov2;
  int l = 0;
  int ngaps;
  int score;
  string gt;
  int gtq;
  bool imprecise;
  // Junction detail for split-read calls, following Manta/DRAGEN conventions.
  // A deletion whose two segments are not flush carries either novel bases
  // between them (ins_len/ins_seq -> SVINSLEN/SVINSSEQ) or an overlap where the
  // shared bases align to both sides (hom_len/hom_seq -> HOMLEN/HOMSEQ). These
  // are reported ALONGSIDE SVLEN/END, never netted out of them: SVLEN stays the
  // number of deleted reference bases and END stays POS + that span.
  int ins_len = 0;
  int hom_len = 0;
  string ins_seq;
  string hom_seq;
  string cigar;
  string reads;
  string sa_reads;
  string rvec;

  SV();
  SV(const string type_, const string &chrom_, uint s_, const string &refall_,
     const string &altall_, const uint w_, const uint cov_, const int ngaps_,
     const int score_, bool imprecise_ = false, uint l_ = 0,
     string cigar_ = ".");
  void add_reads(const vector<string> &reads_);
  // ID of the record holding the other breakend of this junction (BND only).
  // Filled by Caller::link_bnd_mates() once every record exists; stays empty
  // when only one side of the junction cleared the gates, which is the honest
  // state — a dangling MATEID would claim a record that was never written.
  string mate_id;
    void add_sa_reads(const vector<string> &reads_);

  void set_cov(int, int, int, int);
  void set_rvec(const vector<tuple<int, int>> &);
  void set_gt(const string &, int);

  bool operator<(const SV &c) const {
    if (chrom < c.chrom) {
      return true;
    } else if (chrom > c.chrom) {
      return false;
    } else {
      return s < c.s;
    }
  }

  bool operator==(const SV &c) const {
    return chrom == c.chrom and s == c.s and e == c.e and type == c.type;
  }

  friend ostream &operator<<(ostream &os, const SV &sv);
};

#endif
