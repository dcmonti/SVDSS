#include "sfs.hpp"

// TODO: if we split the .sfs into more files, we can load it using more threads
// (but for SVs having a single file shouldn't be the bottleneck)
SFSData parse_sfsfile(const string &sfs_path) {
  spdlog::info("Loading SFSs from {}..", sfs_path);
  SFSData result;
  result.source = SFSSource::UNKNOWN;
  unordered_map<string, vector<SFS>> &SFSs = result.sfss;
  int total = 0;
  string line;
  ifstream inf(sfs_path);
  if (inf.is_open()) {
    string info[4];
    string read_name;
    while (getline(inf, line)) {
      // Parse header lines (starting with #)
      if (!line.empty() && line[0] == '#') {
        if (line.find("#source=bam") != string::npos)
          result.source = SFSSource::BAM;
        else if (line.find("#source=fasta") != string::npos)
          result.source = SFSSource::FASTA;
        continue;
      }
      stringstream ssin(line);
      int i = 0;
      while (ssin.good() && i < 4)
        ssin >> info[i++];
      if (info[0].compare("*") != 0) {
        read_name = info[0];
        SFSs[read_name] = vector<SFS>();
      }
      SFSs[read_name].push_back(
          SFS(read_name, stoi(info[1]), stoi(info[2]), stoi(info[3])));
      ++total;
    }
  }
  if (result.source == SFSSource::UNKNOWN)
    spdlog::warn("SFS file {} has no #source header. Assuming FASTA (legacy behavior). "
                 "Regenerate with latest SVDSS search to add the header.", sfs_path);
  spdlog::info("Loaded {} SFSs from {} reads (source: {}).", total, SFSs.size(),
               result.source == SFSSource::BAM ? "bam" :
               result.source == SFSSource::FASTA ? "fasta" : "unknown (assuming fasta)");
  return result;
}
