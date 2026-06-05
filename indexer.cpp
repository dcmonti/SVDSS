#include "indexer.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <spdlog/spdlog.h>

#include "chromosomes.hpp"

extern "C" {
int main_build(int argc, char *argv[]);
}

using namespace std;

namespace {

// ropebwt3 build short options that take a value (see ropebwt3 build.c usage).
static const string WITH_ARG = "mtplniosS";

bool file_exists(const string &p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0;
}

void touch_file(const string &p) {
  FILE *f = fopen(p.c_str(), "w");
  if (f != nullptr)
    fclose(f);
}

} // namespace

int run_index(int argc, char **argv) {
  bool hpc_mode = false;
  int hpc_cap = 1;
  string output_fmd;
  string append_idx;
  vector<char *> new_argv;
  vector<string> tmp_files;
  vector<char *> owned_strings; // strdup'd entries we must free

  // Pre-scan for HPC settings so positional FASTA compression below uses the
  // right cap regardless of argument order.
  for (int i = 1; i < argc; ++i) {
    string tok = argv[i];
    if (tok == "--hpc")
      hpc_mode = true;
    else if (tok == "--hpc-cap" && i + 1 < argc)
      hpc_cap = atoi(argv[i + 1]);
    else if (tok.rfind("--hpc-cap=", 0) == 0)
      hpc_cap = atoi(tok.c_str() + 10);
  }
  if (hpc_cap < 1)
    hpc_cap = 1;

  new_argv.push_back(argv[0]); // "index"

  for (int i = 1; i < argc; ++i) {
    string tok = argv[i];
    if (tok == "--hpc") {
      continue; // consumed in pre-scan
    }
    if (tok == "--hpc-cap") {
      ++i; // skip the value token too
      continue;
    }
    if (tok.rfind("--hpc-cap=", 0) == 0) {
      continue;
    }
    // option that takes a value (forms: "-X VAL" or "-XVAL")
    if (tok.size() >= 2 && tok[0] == '-' &&
        WITH_ARG.find(tok[1]) != string::npos) {
      string val;
      bool inline_val = tok.size() > 2;
      if (inline_val) {
        val = tok.substr(2);
        new_argv.push_back(argv[i]);
      } else {
        new_argv.push_back(argv[i]);
        if (i + 1 >= argc) {
          spdlog::critical("Missing value for option {}", tok);
          return 1;
        }
        ++i;
        val = argv[i];
        new_argv.push_back(argv[i]);
      }
      if (tok[1] == 'o')
        output_fmd = val;
      else if (tok[1] == 'i')
        append_idx = val;
      continue;
    }
    // bare flag or unknown option -> pass through
    if (tok.size() >= 1 && tok[0] == '-') {
      new_argv.push_back(argv[i]);
      continue;
    }
    // positional FASTA input
    if (hpc_mode) {
      char tmpl[] = "/tmp/svdss_hpc_XXXXXX.fa";
      int fd = mkstemps(tmpl, 3);
      if (fd < 0) {
        spdlog::critical("mkstemps failed for HPC temp file");
        return 1;
      }
      close(fd);
      string tmp_path(tmpl);
      try {
        hpc_write_reference_fasta(string(argv[i]), tmp_path, hpc_cap);
      } catch (const std::exception &e) {
        spdlog::critical("HPC compression of {} failed: {}", argv[i], e.what());
        unlink(tmp_path.c_str());
        return 1;
      }
      tmp_files.push_back(tmp_path);
      char *dup = strdup(tmp_path.c_str());
      owned_strings.push_back(dup);
      new_argv.push_back(dup);
    } else {
      new_argv.push_back(argv[i]);
    }
  }

  if (hpc_mode) {
    if (output_fmd.empty()) {
      spdlog::critical(
          "--hpc requires -o <fmd> so a sidecar marker can be written");
      for (const auto &p : tmp_files)
        unlink(p.c_str());
      for (auto s : owned_strings)
        free(s);
      return 1;
    }
    if (!append_idx.empty() && !file_exists(append_idx + ".hpc")) {
      spdlog::critical(
          "--hpc with -i {}: existing index has no .hpc sidecar; cannot "
          "append non-HPC index in HPC mode",
          append_idx);
      for (const auto &p : tmp_files)
        unlink(p.c_str());
      for (auto s : owned_strings)
        free(s);
      return 1;
    }
    spdlog::info("HPC mode (cap={}): reference inputs homopolymer-compressed "
                 "before indexing",
                 hpc_cap);
  }

  spdlog::info("We rely on 'ropebwt3 build' - just exposing it for convenience");
  int rc = main_build(static_cast<int>(new_argv.size()), new_argv.data());

  for (const auto &p : tmp_files)
    unlink(p.c_str());
  for (auto s : owned_strings)
    free(s);

  if (hpc_mode && rc == 0) {
    touch_file(output_fmd + ".hpc");
    spdlog::info("HPC marker written: {}.hpc", output_fmd);
  }
  return rc;
}
