#ifndef CALLER_HPP
#define CALLER_HPP

#include <ctime>
#include <dirent.h>
#include <iostream>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <abpoa.h>
#include <ksw2.h>
#include <rapidfuzz/fuzz.hpp>
#include <spdlog/spdlog.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "clipper.hpp"
#include "clusterer.hpp"
#include "config.hpp"
#include "sfs.hpp"
#include "sv.hpp"

using namespace std;

// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char _char26_table[256] = {
    0, 1,         2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4 /*'-'*/, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
    4, 1,         4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

class Consensus {
public:
  string seq;
  string chrom;
  string cigar;
  int s;
  int e;

  Consensus(const string _seq, const string _cigar, const string _chrom, int _s,
            int _e) {
    seq = _seq;
    cigar = _cigar;
    chrom = _chrom;
    s = _s;
    e = _e;
  }

  friend ostream &operator<<(ostream &os, const Consensus &c) {
    os << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "\t"
       << "0"
       << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
       << "60"
       << "\t" << c.cigar << "\t"
       << "*"
       << "\t"
       << "0"
       << "\t"
       << "0"
       << "\t" << c.seq << "\t"
       << "*";
    return os;
  }
};

// One POA consensus plus the indices (into the input seqs vector) of the reads
// that abPOA grouped under it. With diploid_poa a heterozygous cluster yields
// two PoaCons (one per allele); otherwise a single one carrying all reads.
struct PoaCons {
  string seq;
  vector<int> read_ids;
};

class Caller {

public:
  void run();

private:
  Configuration *config;
  unordered_map<string, vector<SFS>> SFSs;

  vector<SV> svs;
  vector<Consensus> alignments;

  void pcall(const vector<Cluster> &);
  void clean_dups();
  void filter_sv_chains();
  void write_vcf();
  void write_sam();
  int count_overlapping_original_sfs(const Cluster &cluster,
                                     const Cluster &subcluster,
                                     const SV &sv) const;

  vector<Cluster> split_cluster_by_len(const Cluster &);
  vector<Cluster> split_cluster(const Cluster &);
  vector<PoaCons> run_poa(const vector<string> &);
  // Align query to ref with the same ksw2 scoring/mode as the tumor consensus;
  // returns the CIGAR as (length, op) pairs (op: 0=M, 1=I, 2=D).
  vector<pair<int, int>> ksw_align(const string &ref, const string &query);
  // Re-align a normal read's window around sv with ksw_align and test for a
  // same-type indel of concordant size (used by --germline-realign).
  bool normal_read_concordant(bam1_t *aln, const string &chrom, const SV &sv,
                              int cl_s, int cl_e, int want, float min_ratio_len,
                              int diff_max, int diff_min_len);
  // Count primary tumour reads spanning [s, s+len] with a single CIGAR D op
  // reciprocally overlapping the interval (>= 0.5). These "through-reads" are
  // minimap2's non-clipped representation of a < ~15 kbp deletion and are
  // read-disjoint from the split (clipped) reads (Mode B clip rescue).
  uint count_through_del_reads(samFile *bam, hts_idx_t *idx, bam_hdr_t *hdr,
                               const string &chrom, uint s, uint len);

  // parallelize
  vector<vector<SV>> _p_svs;
  vector<vector<Consensus>> _p_alignments;
  // per-thread germline-filtered regions: {sv.chrom, sv.s, sv.s + abs(sv.l)}
  // populated in pcall(); merged into the per-chrom clip vartrees in run()
  vector<vector<std::tuple<std::string, int, int>>> _p_germline_regions;

  // per-thread handles for normal contigs BAM (optional germline filter)
  samFile **_p_normal_bam = nullptr;
  hts_idx_t **_p_normal_idx = nullptr;
  bam_hdr_t **_p_normal_hdr = nullptr;

  bool is_germline(const SV &sv, const string &chrom, int cl_s, int cl_e,
                   int t);
  // Germline check for breakend-like events (BND/INV) using the split (SA)
  // alignments of the normal contigs instead of CIGAR I/D ops.
  bool is_germline_breakend(const SV &sv, int t);

  // True if the reference carries an N within +-window bp of pos on chrom
  // (0-based access into chromosome_seqs; pos is the 1-based VCF position).
  bool ref_has_N_near(const string &chrom, long pos, int window) const;
  // Parse a BND ALT ("N]chr5:1234]") into mate chrom/pos. False if not a BND.
  bool bnd_mate(const SV &sv, string &mate_chrom, long &mate_pos) const;
  // Combined exclusion: bed overlap or N-proximity at the primary breakpoint
  // (and, for BND, at the mate breakpoint too).
  bool excluded_by_bed_or_N(const SV &sv) const;

  void print_vcf_header();
};

#endif
