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
  // Only set_gt() ever writes gtq, and only the POA path calls it. The clipped
  // path builds its records straight from the constructor, so without a default
  // here every clipped call printed whatever happened to be on the stack: on
  // HG008 that read 0 for 108 of 116 clipped records and garbage (1,
  // -2010473030, 1886501328) for the rest, which also made two runs of the same
  // binary disagree. 0 is the honest value for a path that computes no
  // genotype quality, and matches the neutral gt/rvec the constructor sets.
  int gtq = 0;
  bool imprecise;
  // Junction detail for split-read calls, following Manta/DRAGEN conventions.
  // A deletion whose two segments are not flush carries either novel bases
  // between them (ins_len/ins_seq -> SVINSLEN/SVINSSEQ) or an overlap where the
  // shared bases align to both sides (hom_len/hom_seq -> HOMLEN/HOMSEQ). These
  // are reported ALONGSIDE SVLEN/END, never netted out of them: SVLEN stays the
  // number of deleted reference bases and END stays POS + that span.
  int ins_len = 0;
  int hom_len = 0;
  // Reference bases replaced by the insertion (clipped path only). A clipped INS
  // is a replacement event: dQ query bases stand in for dR reference bases, and
  // SVLEN carries only the net dQ - dR. Without dR the partner breakpoint is not
  // recoverable from the record, and the germline filter — which has to know
  // where the mate segment resumes — would assume it sits at POS. 0 for a pure
  // insertion (POA path, or a flush clipped junction).
  int ref_gap = 0;
  // Deletions that the consensus alignment split into several D operations and
  // the caller put back together. del_parts counts them (0 or 1 = untouched
  // record), del_kept is how many reference bases between them the alignment
  // does NOT delete. Additive, like the fields above: SVLEN stays the full
  // POS..END span, so nothing downstream has to learn about the split, and
  // del_kept says how much of that span the alignment actually kept. The exact
  // structure remains readable in CIGAR.
  int del_parts = 0;
  int del_kept = 0;
  string ins_seq;
  string hom_seq;
  string cigar;
  string reads;
  string sa_reads;
  string rvec;
  // ID of the record holding the other breakend of this junction (BND only).
  // Filled by Caller::link_bnd_mates() once every record exists; stays empty
  // when only one side of the junction cleared the gates, which is the honest
  // state — a dangling MATEID would claim a record that was never written.
  string mate_id;

  SV();
  SV(const string type_, const string &chrom_, uint s_, const string &refall_,
     const string &altall_, const uint w_, const uint cov_, const int ngaps_,
     const int score_, bool imprecise_ = false, uint l_ = 0,
     string cigar_ = ".");
  void add_reads(const vector<string> &reads_);
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
