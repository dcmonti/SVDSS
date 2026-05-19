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
                     int *hpc_to_orig) {
  if (len <= 0) {
    hpc_to_orig[0] = 0;
    return 0;
  }
  int j = 0;
  dst[0] = src[0];
  hpc_to_orig[0] = 0;
  for (int i = 1; i < len; ++i) {
    if (src[i] != src[i - 1]) {
      ++j;
      dst[j] = src[i];
      hpc_to_orig[j] = i;
    }
  }
  int hpc_len = j + 1;
  hpc_to_orig[hpc_len] = len; // sentinel: one past the end of the last run
  return hpc_len;
}

void hpc_write_reference_fasta(const string &in_fa, const string &out_fa) {
  spdlog::info("HPC-compressing reference {} -> {}..", in_fa, out_fa);
  gzFile fp = gzopen(in_fa.c_str(), "r");
  if (fp == nullptr)
    throw runtime_error("Cannot open reference " + in_fa);
  kseq_t *seq = kseq_init(fp);
  FILE *out = fopen(out_fa.c_str(), "w");
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
    // HPC-compress in place, uppercased
    char prev = '\0';
    int col = 0;
    int kept = 0;
    for (int i = 0; i < l; ++i) {
      char c = (char)toupper((unsigned char)seq->seq.s[i]);
      if (c == prev)
        continue;
      prev = c;
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
