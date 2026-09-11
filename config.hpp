#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <string>
#include <vector>

#include "cxxopts.hpp"
#include "bed_filter.hpp"

using namespace std;

// clang-format off

static const char MAIN_USAGE_MESSAGE[] =
  "SVDSS [index|smooth|search|call] --help\n"
  "      --help                           print help message\n"
  "      --version                        print version\n";

static const char INDEX_USAGE_MESSAGE[] =
  "SVDSS index --reference <reference> --index <index>\n"
  "      --binary                         store the index in binary format\n"
  "      --append <oldindex>              append to existing index\n"
  "      --help                           print help message\n";

static const char SMOOTH_USAGE_MESSAGE[] =
  "SVDSS smooth --reference <reference> --bam <bam>\n"
  "      --min-mapq <INT>                 minimum mapping quality (default: 20)\n"
  "      --min-indel-length <INT>         indels shorter than this are smoothed out (default: 20)\n"
  "      --accp <FLOAT>                   accuracy percentile for alignments filtering (default: 0.98)\n"
  "      --germline-min-ro <FLOAT>        reciprocal-overlap to call an SV germline vs normal contigs (default: 0.95)\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

static const char SEARCH_USAGE_MESSAGE[] =
  "SVDSS search --index <index> [--bam <bam>|--fastq <fastq>]\n"
  "      --bsize <INT>                    batch size (default: 10000)\n"
  "      --noputative                     when input is smoothed bam, do not filter unsmoothed reads (default: putative)\n"
  "      --noassemble                     do not assemble specific strings overlapping on a read (default: assemble)\n"
  "      --hpc                            homopolymer-compress reads before search; the index must have been built with --hpc\n"
  "      --hpc-cap <INT>                  cap each homopolymer run at INT copies (default: 1 = full collapse); must match the index\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

static const char CALL_USAGE_MESSAGE[] =
  "SVDSS call --reference <reference> --bam <bam> --sfs <sfs>\n"
  "      --poa <FILE>                     store POA in .sam format to this file (default: do not store)\n"
  "      --clusters <FILE>                store clusters to this file (default: do not store)\n"
  "      --clip-clusters <FILE>           store clipped-SFS clusters to this file (requires --clipped; default: do not store)\n"
  "      --min-cluster-weight <INT>       minimum number of supporting superstrings for a call to be reported (default: 2)\n"
  "      --min-sv-length <INT>            minimum length of reported SVs (default: 25)\n"
  "      --require-sfs-overlap            require at least one original SFS to overlap the SV interval\n"
  "      --noht                           do not use haplotagging information even if present\n"
  // "      --noref                          do not report 0/0 calls\n"
  "      --min-mapq                       minimum mapping quality (default: 20)\n"
  "      --max-cluster-dist <INT>          maximum distance for clustering SFSs (default: 0, auto-computed)\n"
  "      --clipped                        calls SVs from clipped SFS (EXPERIMENTAL)\n"
  "      --no-clipped-fallback            disable paired-clip fallback calls (default: false)\n"
  "      --normal-contigs-bam <FILE>      BAM of normal-tissue contigs aligned to reference for germline filtering\n"
  "      --germline-min-ro <FLOAT>        length-ratio to call an SV germline vs normal reads (default: 0.95)\n"
  "      --germline-diff-max <INT>        also call germline if |sv_len - normal_indel_len| < this (0=min_sv_length, default: 0)\n"
  "      --germline-diff-min-len <INT>    minimum normal-indel length for the difference-based germline test (default: 10)\n"
  "      --no-germline-diff               disable the difference-based germline test (ratio-only)\n"
  "      --diploid-poa                    emit up to two POA consensuses per cluster (one per allele; experimental)\n"
  "      --poa-min-freq <FLOAT>           min allele fraction for --diploid-poa to split a cluster (default: 0.20)\n"
  "      --no-gates                       disable the post-call GERMLINE/LOWPOWER gates\n"
  "      --gate-min-reads <INT>           normal reads matching ALEN needed to drop a call as germline (default: 2)\n"
  "      --gate-tol-frac <FLOAT>          relative tolerance of the ALEN match (default: 0.25)\n"
  "      --gate-q <FLOAT>                 P(a normal read of the event's haplotype shows a het) (default: 0.85)\n"
  "      --gate-conf <FLOAT>              below this confidence a call is flagged LOWPOWER (default: 0.95)\n"
  "      --germline-realign               re-align normal reads with the tumor's ksw2 for the germline test (experimental)\n"
  "      --germline-realign-margin <INT>  flank added to the SV window when re-aligning normal reads (default: 300)\n"
  "      --no-merge-del                   do not merge D operations of the same alignment split by a short aligned stretch\n"
  "      --merge-del-gap <INT>            maximum gap between merged D operations, in bp (default: 100)\n"
  "      --merge-del-gap-frac <FLOAT>     maximum gap as a fraction of the deleted bases so far (default: 0.20)\n"
  "      --no-merge-del-xc                do not rejoin deletions split across consensus alignments\n"
  "      --merge-del-xc-gap <INT>         maximum gap for that rejoin, in bp (default: 5000)\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

// clang-format on

class Configuration {

private:
  static Configuration *instance;

public:
  static Configuration *getInstance();

  void parse(int argc, char *argv[]);
  void print_help(const string &) const;

  // general
  int threads = 4;
  int batch_size = 10000;
  bool version = false;
  bool verbose = false;
  bool help = false;

  // smoother
  float accp = 0.98;

  // pingpong.index
  bool binary = false;
  // pingpong.search
  bool assemble = true;
  bool putative = true;
  bool hpc = false; // homopolymer-compressed search (requires HPC-built index)
  int hpc_cap = 1;  // cap per homopolymer run (1 = full collapse); must match index
  int overlap = -1;
  int max_output = 100000;
  // call
  uint flank = 100;
  uint ksize = 7;
  uint min_sv_length = 25;
  uint min_mapq = 20;
  int min_indel_length = 20;
  uint min_cluster_weight = 2;
  // Clipped/SA calling: same-chromosome, OPPOSITE-strand split alignments whose
  // two ends are at least this far apart (bp) are treated as an intra-chromosomal
  // translocation and emitted as a BND breakend (like cross-chrom), instead of one
  // giant INV spanning the whole gap. Same-strand events (DEL/DUP/INS) are never
  // affected. 0 disables (keep legacy behaviour).
  uint min_bnd_dist = 5000000;
  float min_ratio = 0.97; // FIXME: change name
  // Length-concordance threshold for the germline filter: an SV is called
  // germline if normal reads (or contigs) carry a matching event whose length
  // ratio with the SV is >= this. Lower it to filter more germline (e.g. VNTR
  // length-polymorphisms).
  float germline_min_ro = 0.90;
  // Read-based germline filter: minimum number of concordant normal reads to
  // declare an SV germline, and the maximum number of normal reads examined per
  // SV (bounds cost in high-depth/repeat loci). Reads below min_mapq are
  // skipped. The filter queries config.normal_contigs_bam, which may point at a
  // normal *reads* BAM instead of contigs.
  // min_reads=1: one concordant normal read is enough, because a somatic event
  // has no normal support at all -- a single concordant read contradicts the
  // somatic hypothesis rather than weakly disfavouring it.
  int germline_min_reads = 1;
  int germline_max_reads = 1000;
  // Difference-based germline concordance (union with the ratio test above). A
  // normal read also counts as concordant when it carries a same-type indel of
  // length L_n >= germline_diff_min_len with |sv_len - L_n| < germline_diff_max.
  // This catches VNTR/low-complexity length-jitter that the strict ratio test
  // misses (normal allele present but a different number of repeat units), while
  // leaving true somatic events untouched (no normal support at all). Enabled by
  // default; --no-germline-diff restores the ratio-only behaviour. A
  // germline_diff_max of 0 means "use min_sv_length" (resolved in Caller).
  bool germline_diff = true;
  int germline_diff_max = 0; // 0 -> min_sv_length
  int germline_diff_min_len = 10;
  // Diploid POA: let abPOA emit up to two consensus sequences (one per allele)
  // per cluster, so heterozygous/multi-allele VNTR loci are not averaged into a
  // phantom-size consensus. poa_min_freq is abPOA's min allele fraction to split.
  // On by default; --no-diploid-poa restores the single-consensus behaviour.
  bool diploid_poa = true;
  float poa_min_freq = 0.20;
  // Germline re-alignment: for the germline test, align each normal read's local
  // window with the SAME ksw2 as the tumor consensus instead of reading the BAM's
  // minimap2 CIGAR. In repeats the two aligners place gaps differently (ksw2
  // fragments a deletion, minimap2 keeps it whole), so a tumor fragment finds no
  // match among the BAM CIGARs and evades the filter. Re-aligning makes both
  // sides fragment identically. On by default (--no-germline-realign disables).
  // Additive: it can only mark MORE reads concordant, so it can only remove
  // calls. Skipped for SVs longer than germline_realign_max_len (window too
  // large to realign cheaply).
  // Post-call gates, evaluated once on the records that are about to be written
  // (see Caller::apply_gates). They do NOT replace the germline filter that runs
  // during calling.
  //   GERMLINE: >= gate_min_reads normal reads carry the same allele -> kept,
  //             FILTER=GERMLINE. FLAGGED, not dropped: the threshold is
  //             empirical, so the record stays auditable and
  //             `bcftools view -f PASS` gives the same set a hard drop would.
  //   LOWPOWER: too little normal depth ON THE EVENT'S HAPLOTYPE for "absent
  //             from the normal" to mean anything -> kept, FILTER=LOWPOWER.
  bool gates = true;
  int gate_min_reads = 2;
  float gate_tol_frac = 0.25;
  // gate_win: half-window, in bp, over which gate A looks for normal reads
  //   carrying a same-type indel of MATCHING length (it also sets the fetch
  //   region and the alen_med search). Was hardcoded at 150 in
  //   collect_call_stats; measured too tight on both samples. On COLO829 the
  //   germline evidence sits 152, 185 and 368 bp from the breakpoint, and at
  //   chr11:10249845 82 of 84 spanning normal reads carry the allele while the
  //   gate at 150 saw none -- missed by 2 bp. Widening to 500 marked 5 more
  //   COLO829 records and 4 more HG008 records GERMLINE, ALL of them false
  //   positives, and cost 0 true positives on either sample.
  //   Deliberately a SEPARATE knob from gate_poly_win below: the two gates ask
  //   different questions and must be tunable apart, even though both currently
  //   sit at 500. The fetch region uses max(gate_win, gate_poly_win).
  int gate_win = 500;
  // ---- gate C: locus length-polymorphic in the normal ----------------------
  // Gate A asks whether the normal carries THIS allele, so it compares lengths.
  // At a length-polymorphic locus that test cannot fire: the normal allele
  // differs in length by definition. Gate C drops the length comparison and asks
  // instead what FRACTION of spanning normal reads carries a same-type indel at
  // all -- "is this locus polymorphic in the normal", which is the question that
  // actually decides whether a tumour call of some other length is somatic.
  //
  // gate_poly_win: half-window, in bp, for collecting those indels. Unlike
  //   gate_win this window is CENTRED ON THE BREAKPOINT and does not span the
  //   event: with the spanning window, truthset_22 -- a 32.5 kb COLO829
  //   deletion -- was flagged POLYNORMAL because within 32 kb the normal has
  //   unrelated germline indels (23 reads of 64), costing a true positive.
  // gate_poly_frac: fraction above which the locus is called polymorphic.
  //   0.40, not 0.25, and the reason matters. Calibrating against truvari labels
  //   alone put the plateau at 0.15-0.30, but those labels call two HG008
  //   records false that manual screening keeps: INS chr15:74894107 (0.333, a
  //   174 bp expansion of a ~90 bp germline allele) and INS chr4:30278703
  //   (0.364, a 28 bp-motif VNTR contraction). Both are in the curated
  //   keep-list, so the GT and the intended behaviour disagree exactly there.
  //   There is a clean gap: those two sit at 0.333 and 0.364, while every
  //   COLO829 false positive sits at 0.4375, 0.4545 and 1.00. Any threshold in
  //   (0.364, 0.4375] keeps both and still catches all three; 0.40 is near the
  //   middle of that window.
  //   Measured at 0.40 on TWO samples. COLO829: all 3 POLYNORMAL still caught,
  //   F1 unchanged. HG008: 9 of the FP caught, and the one true positive given
  //   up is chr1:224013772 (0.565, a 2490 bp insertion where 13 of 23 normal
  //   reads already carry one -- no threshold saves it without giving up most
  //   of the FP). RE-MEASURE on a third sample.
  // gate_poly_min_reads: absolute floor. A fraction over 4 reads is noise: some
  //   HG008 loci gave 2/4 = 0.50. Higher than gate_min_reads on purpose, because
  //   gate A's length-matched test is far more specific than this one.
  int gate_poly_win = 500;
  float gate_poly_frac = 0.40;
  int gate_poly_min_reads = 3;
  // ---- gate D: mismapping sink --------------------------------------------
  // Gate C is blind by dilution where the denominator is absurd: at
  // chr16:46391338 the normal has 25568 spanning reads of which 263 carry the
  // indel, a fraction of 0.01 -- yet 263 reads of germline evidence is a
  // mountain. True positives on COLO829 span 16-98 normal reads, so 25568 is
  // ~300x the norm: not a locus, a satellite mismapping sink.
  //
  // Self-normalising on purpose: the MEDIAN of n_exam over the records of the
  // run estimates this sample's normal depth, so there is no extra BAM pass and
  // no external annotation, and it follows the coverage. Between 98 and 23910
  // any factor from 2 to 200 separates, so 10 is deliberately loose. It reads
  // the NORMAL's depth, so tumour copy-number gain cannot trigger it.
  int gate_depth_factor = 10;
  float gate_q = 0.85;   // P(a normal read of the right haplotype shows a het)
  float gate_conf = 0.95;
  bool germline_realign = true;
  int germline_realign_margin = 300; // (unused: window is now the cluster interval)
  int germline_realign_max_len = 15000; // skip realign when cl_e-cl_s exceeds this
  bool useht = true;
  // bool noref = false;
  bool clipped = false;
  bool clipped_fallback = true;
  // Require at least 2 of the original (pre-assembly) SFS to overlap the SV
  // interval. Applies to the SFS/POA path only — clipped-path SVs carry no SFS
  // and are never touched. On by default; --no-require-sfs-overlap disables.
  bool require_sfs_overlap = true;
  // Merge D operations of the SAME consensus alignment that a short aligned
  // stretch splits apart. ksw2 fragments a deletion where minimap2 keeps it
  // whole (the same effect germline_realign compensates for on the normal side),
  // so one event leaves the caller as two records that each miss the truth on
  // size. Scoped to the POA path and to DEL: clipped-path calls carry no CIGAR
  // and never enter, and for INS the length lives in the sequence rather than in
  // the reference span, so a gap-based rule there is meaningless.
  // Both gates must hold, the relative one keeping small neighbouring deletions
  // apart: two 60 bp DELs 100 bp away fail it (100 > 0.2 * 120).
  bool merge_del = true;
  int merge_del_max_gap = 100;
  float merge_del_max_gap_frac = 0.2;
  // Cross-consensus rejoin, keyed on the supporting reads being the SAME set.
  // The gap cap is a backstop only: read-set identity is what makes it safe.
  bool merge_del_xc = true;
  int merge_del_xc_max_gap = 5000;
  int max_cluster_dist = 0; // 0 means auto-computed

  string bam = "";
  string sfs = "";
  string append = "";
  string index = "";
  string reference = "";
  string fastq = "";
  string poa = "";
  string clusters = "";
  string clip_clusters = "";
  string normal_contigs_bam = "";
  string bed_exclusion = "";
  
  BedFilter bed_filter;

private:
  Configuration();

  Configuration(Configuration const &) = delete;
  void operator=(Configuration const &) = delete;

  Configuration &operator[](string);

  cxxopts::Options parser;
};

#endif
