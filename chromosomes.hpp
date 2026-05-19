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
// seq_nt6_table in ping_pong.hpp). Collapses each maximal run of identical
// codes to a single code. Writes the compressed sequence to `dst` (which must
// have capacity >= len) and fills `hpc_to_orig` (length >= len + 1) with the
// original positions of each compressed run, plus a sentinel at index hpc_len
// equal to `len`. Returns the compressed length.
//
// Example (codes shown as chars for clarity):
//   src         = AAAATTGCC      (len=9)
//   dst         = ATGC           (hpc_len=4)
//   hpc_to_orig = [0, 4, 6, 7, 9]
// So an HPC interval [begin=1, end=2] (covering T,G) translates back to
// original [hpc_to_orig[1]=4, hpc_to_orig[3]-1=6] = TTG.
int hpc_compress_nt6(const uint8_t *src, int len, uint8_t *dst,
                     int *hpc_to_orig);

// Read a FASTA from `in_fa` and write the homopolymer-compressed version to
// `out_fa`. Used by the index step to feed ropebwt3 a compressed reference.
// Operates on plain ASCII (A,C,G,T,N) preserving case-folded uppercase.
void hpc_write_reference_fasta(const string &in_fa, const string &out_fa);

#endif
