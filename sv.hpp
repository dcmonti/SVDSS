#ifndef SV_HPP
#define SV_HPP

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Per-record evidence for telling a somatic call from a germline one, emitted
// as INFO so the decision can be prototyped on the VCF before any of it becomes
// a filter. Nothing here changes which records are written.
//
// Two independent axes, deliberately kept apart:
//
//   1. Is the normal carrying this event?  n_conc/n_lens say what the normal
//      actually holds. is_germline() answers only yes/no against sv.l, so a
//      record whose reported length differs from the allele the reads carry
//      escapes with no trace of why; n_lens keeps the raw lengths so the test
//      can be re-run against any length downstream.
//
//   2. Could we even have seen it?  A germline het sits on ONE haplotype, so
//      the power to exclude it is set by the normal depth on THAT haplotype,
//      not by the total. nhp1/nhp2 vs hp1/hp2 make that computable: a locus can
//      carry several normal reads all from the haplotype opposite the event,
//      where total depth alone looks like adequate evidence of absence.
//
// -1 means "not evaluated" (no normal BAM given) and must not be read as zero.
// n_lens/n_conc stay empty for BND/DUP/INV, which carry no indel to compare;
// the haplotype counts are filled for every type.
struct CallStats {
  // Tumour side: haplotypes of the reads supporting THIS record, and the indel
  // length those reads actually carry. A record built from reads of both
  // haplotypes is not a single-copy event; a wide alen_min-alen_max means the
  // consensus is averaging reads that do not agree on a length.
  //
  // alen is measured the same way NLENS is measured on the normal, so the two
  // are directly comparable -- which SVLEN and NLENS are not, because SVLEN
  // comes out of the consensus alignment and can differ from every read that
  // built it: a DEL can be reported at one length while its reads carry another.
  int hp1 = 0, hp2 = 0, hp0 = 0;
  int alen_med = 0, alen_min = 0, alen_max = 0;
  int n_sup = 0; // supporting reads actually found in the tumour BAM
  // Normal side, filled by Caller::collect_call_stats().
  int n_exam = -1;    // reads passing flag/MAPQ and spanning the breakpoint
  int n_lowq = -1;    // reads dropped for MAPQ alone (min_mapq is strict)
  int n_conc = -1;    // reads carrying a concordant indel (germline support)
  int nhp1 = -1, nhp2 = -1, nhp0 = -1; // spanning normal reads per haplotype
  vector<int> n_lens; // same-type indel lengths seen in the normal window
};

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
  // here every clipped call printed whatever happened to be on the stack, which
  // also made two runs of the same binary disagree. 0 is the honest value for a
  // path that computes no
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
  CallStats stats;
  // VCF FILTER. "PASS" unless a post-call gate demoted it: a flagged record is
  // still written, so the full file stays the ALL callset and `bcftools view -f
  // PASS` is the confident subset.
  string filter = "PASS";

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
