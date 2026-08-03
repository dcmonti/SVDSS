#include "chromosomes.hpp"

#include <cstdio>
#include <cstring>
#include <stdexcept>

KSEQ_INIT(gzFile, gzread)

// FIXME: avoid this
vector<string> chromosomes;
unordered_map<string, char *> chromosome_seqs;

void load_chromosomes(string path) {
  spdlog::info("Loading reference genome from {}..", path);
  gzFile fp = gzopen(path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    string name(seq->name.s);
    spdlog::debug("Extracted {} ({}bp)", seq->name.s, l);
    for (int i = 0; i < l; i++)
      seq->seq.s[i] = toupper(seq->seq.s[i]);
    chromosomes.push_back(seq->name.s);
    char *s = (char *)malloc(sizeof(char) * (l + 1));
    memcpy(s, seq->seq.s, l + 1);
    s[l] = '\0';
    chromosome_seqs[seq->name.s] = s;
  }
  kseq_destroy(seq);
  gzclose(fp);
}

void destroy_chromosomes() {
  for(const auto &chrom : chromosomes) {
    free(chromosome_seqs[chrom]);
  }
}

int hpc_compress_nt6(const uint8_t *src, int len, uint8_t *dst,
                     int *hpc_to_orig, int cap) {
  if (cap < 1)
    cap = 1;
  if (len <= 0) {
    hpc_to_orig[0] = 0;
    return 0;
  }
  int j = 0;
  int i = 0;
  while (i < len) {
    int run_start = i;
    uint8_t c = src[i];
    int run_end = i + 1;
    while (run_end < len && src[run_end] == c)
      ++run_end;
    int run_len = run_end - run_start;
    int keep = run_len < cap ? run_len : cap;
    // Emit `keep` copies. The first keep-1 each own one original position; the
    // last one owns the remaining tail of the run, so the slices stay
    // contiguous and the sentinel below closes the partition.
    for (int t = 0; t < keep; ++t) {
      dst[j] = c;
      hpc_to_orig[j] = run_start + t;
      ++j;
    }
    i = run_end;
  }
  int hpc_len = j;
  hpc_to_orig[hpc_len] = len; // sentinel: one past the end of the last run
  return hpc_len;
}

void hpc_write_reference_fasta(const string &in_fa, const string &out_fa,
                               int cap, bool append) {
  if (cap < 1)
    cap = 1;
  spdlog::info("HPC-compressing reference {} -> {} (cap={})..", in_fa, out_fa,
               cap);
  gzFile fp = gzopen(in_fa.c_str(), "r");
  if (fp == nullptr)
    throw runtime_error("Cannot open reference " + in_fa);
  kseq_t *seq = kseq_init(fp);
  FILE *out = fopen(out_fa.c_str(), append ? "a" : "w");
  if (out == nullptr) {
    kseq_destroy(seq);
    gzclose(fp);
    throw runtime_error("Cannot open output " + out_fa);
  }
  static constexpr int LINE_W = 80;
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    // Write FASTA header (with optional comment, mimicking kseq behaviour)
    if (seq->comment.l)
      fprintf(out, ">%s %s\n", seq->name.s, seq->comment.s);
    else
      fprintf(out, ">%s\n", seq->name.s);
    // HPC-compress in place, uppercased: keep at most `cap` copies per run.
    char prev = '\0';
    int run = 0; // copies of the current run already emitted
    int col = 0;
    int kept = 0;
    for (int i = 0; i < l; ++i) {
      char c = (char)toupper((unsigned char)seq->seq.s[i]);
      if (c == prev) {
        if (run >= cap)
          continue; // run already at cap, drop this base
        ++run;
      } else {
        prev = c;
        run = 1;
      }
      fputc(c, out);
      ++kept;
      if (++col == LINE_W) {
        fputc('\n', out);
        col = 0;
      }
    }
    if (col != 0)
      fputc('\n', out);
    spdlog::debug("HPC {}: {}bp -> {}bp ({:.1f}%)", seq->name.s, l, kept,
                  l > 0 ? 100.0 * kept / l : 0.0);
  }
  fclose(out);
  kseq_destroy(seq);
  gzclose(fp);
}
