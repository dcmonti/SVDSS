#ifndef INDEXER_HPP
#define INDEXER_HPP

// Wrapper around ropebwt3's `main_build` that adds an optional
// homopolymer-compressed indexing mode. When `--hpc` is present in argv:
//   - each positional FASTA input is HPC-compressed to a temporary file
//     before being passed to ropebwt3 (the temp file replaces the original
//     path in argv);
//   - the `-o <fmd>` output path is captured; on a successful build, an empty
//     sidecar file `<fmd>.hpc` is written so the search step can detect that
//     the index is HPC and compress reads accordingly;
//   - if `-i <existing>` is also given, `<existing>.hpc` must already exist
//     (cannot mix HPC and non-HPC indices in an append).
//
// All other arguments are forwarded verbatim to ropebwt3's `main_build`.
// argv[0] is expected to be the subcommand name ("index"), matching how
// ropebwt3 ketopt-parses its own argv.
int run_index(int argc, char **argv);

#endif
