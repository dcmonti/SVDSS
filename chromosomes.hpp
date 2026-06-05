#ifndef CHR_HPP
#define CHR_HPP

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include <spdlog/spdlog.h>

#include "config.hpp"
#include "kseq.h"

using namespace std;

extern vector<string> chromosomes;
extern unordered_map<string, char *> chromosome_seqs;

void load_chromosomes(string path);
void destroy_chromosomes();

// Homopolymer-compress an nt6-encoded sequence (alphabet 0..5 from
// seq_nt6_table in ping_pong.hpp). Each maximal run of identical codes is
// capped at `cap` copies (runs shorter than `cap` are left untouched); cap=1
// is the classic full homopolymer compression (run -> single code). Writes the
// compressed sequence to `dst` (capacity >= len) and fills `hpc_to_orig`
// (length >= len + 1) so that each compressed position owns a contiguous slice
// of original positions, with a sentinel at index hpc_len equal to `len`.
// Returns the compressed length. The slices partition [0, len) exactly, so an
// HPC interval [begin, end] maps to original [hpc_to_orig[begin],
// hpc_to_orig[end + 1] - 1].
//
// Example with cap=1 (codes shown as chars for clarity):
//   src         = AAAATTGCC      (len=9)
//   dst         = ATGC           (hpc_len=4)
//   hpc_to_orig = [0, 4, 6, 7, 9]
// Example with cap=2:
//   src         = AAAATTGCC      (len=9)
//   dst         = AATTGCC        (hpc_len=7)
//   hpc_to_orig = [0, 1, 4, 5, 6, 7, 8, 9]
int hpc_compress_nt6(const uint8_t *src, int len, uint8_t *dst,
                     int *hpc_to_orig, int cap = 1);

// Read a FASTA from `in_fa` and write the homopolymer-compressed version to
// `out_fa`, capping each run at `cap` copies (cap=1 = full collapse). Used by
// the index step to feed ropebwt3 a compressed reference. Operates on plain
// ASCII (A,C,G,T,N) preserving case-folded uppercase. The cap must match the
// one used to compress reads at search time.
void hpc_write_reference_fasta(const string &in_fa, const string &out_fa,
                               int cap = 1);

#endif
