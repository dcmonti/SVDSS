#include "config.hpp"

Configuration *Configuration::instance = nullptr;

Configuration *Configuration::getInstance() {
  if (instance == nullptr) {
    instance = new Configuration();
  }
  return instance;
}

void Configuration::print_help(const string &mode) const {
  if (mode.compare("index") == 0) {
    cerr << INDEX_USAGE_MESSAGE << endl;
  } else if (mode.compare("smooth") == 0) {
    cerr << SMOOTH_USAGE_MESSAGE << endl;
  } else if (mode.compare("search") == 0) {
    cerr << SEARCH_USAGE_MESSAGE << endl;
  } else if (mode.compare("call") == 0) {
    cerr << CALL_USAGE_MESSAGE << endl;
  } else {
    cerr << MAIN_USAGE_MESSAGE << endl;
  }
}

Configuration::Configuration()
    : parser(
          "SVDSS, Structural Variant Discovery from Sample-specific Strings") {
  // clang-format off
  parser.add_options()
    ("bam", "", cxxopts::value<std::string>())
    ("sfs", "", cxxopts::value<std::string>())
    ("poa", "", cxxopts::value<std::string>())
    ("clusters", "", cxxopts::value<std::string>())
    ("clip-clusters", "", cxxopts::value<std::string>())
    ("normal-contigs-bam", "", cxxopts::value<std::string>())
    ("index", "", cxxopts::value<std::string>())
    ("fastx", "", cxxopts::value<std::string>())
    ("reference", "", cxxopts::value<std::string>())
    ("append", "", cxxopts::value<std::string>())
    ("threads", "", cxxopts::value<int>())
    ("bsize", "", cxxopts::value<int>())
    ("omax", "", cxxopts::value<int>())
    ("min-sv-length", "", cxxopts::value<int>())
    ("min-indel-length", "", cxxopts::value<int>()) // threshold for smoothing minimum length
    ("min-mapq", "", cxxopts::value<int>())
    ("min-cluster-weight", "", cxxopts::value<int>())
    ("min-bnd-dist", "", cxxopts::value<int>()) // same-chrom opposite-strand split >= this bp -> BND (0 disables)
    ("max-cluster-dist", "", cxxopts::value<int>())
    ("accp", "", cxxopts::value<float>())
    ("germline-min-ro", "", cxxopts::value<float>())
    ("germline-min-reads", "", cxxopts::value<int>()) // read-based germline: min concordant normal reads
    ("germline-max-reads", "", cxxopts::value<int>()) // read-based germline: max normal reads examined per SV
    ("germline-diff-max", "", cxxopts::value<int>()) // difference-based germline: |sv_len-normal_len| < this (0=min_sv_length)
    ("germline-diff-min-len", "", cxxopts::value<int>()) // difference-based germline: min normal-indel length
    ("no-germline-diff", "", cxxopts::value<bool>()->default_value("false")) // disable difference-based germline test
    ("diploid-poa", "", cxxopts::value<bool>()->default_value("false")) // emit up to two POA consensuses per cluster (now ON by default: kept as a no-op for backwards compatibility)
    ("no-diploid-poa", "", cxxopts::value<bool>()->default_value("false")) // single POA consensus per cluster
    ("poa-min-freq", "", cxxopts::value<float>()) // min allele fraction for diploid-poa split
    ("germline-realign", "", cxxopts::value<bool>()->default_value("false")) // re-align normal reads with tumor ksw2 for germline test (now ON by default: kept as a no-op for backwards compatibility)
    ("no-germline-realign", "", cxxopts::value<bool>()->default_value("false")) // germline test on the BAM CIGAR only, no ksw2 fallback
    ("germline-realign-margin", "", cxxopts::value<int>()) // flank added to SV window when re-aligning
    ("germline-realign-max-len", "", cxxopts::value<int>()) // skip realign for SVs larger than this
    ("require-sfs-overlap", "", cxxopts::value<bool>()->default_value("false")) // require >=2 original SFS overlapping the SV interval (now ON by default: kept as a no-op for backwards compatibility)
    ("no-require-sfs-overlap", "", cxxopts::value<bool>()->default_value("false")) // do not require original SFS to overlap the SV interval
    ("bed-exclusion", "", cxxopts::value<std::string>()) // BED file of regions to exclude from analysis. Any SVs overlapping these regions will be filtered out.
    ("clipped", "", cxxopts::value<bool>()->default_value("false"))
    ("no-clipped-fallback", "", cxxopts::value<bool>()->default_value("false"))
    // ("noref", "", cxxopts::value<bool>()->default_value("false"))
    ("noht", "", cxxopts::value<bool>()->default_value("false"))
    ("noassemble", "", cxxopts::value<bool>()->default_value("false"))
    ("noputative", "", cxxopts::value<bool>()->default_value("false"))
    ("hpc", "", cxxopts::value<bool>()->default_value("false"))
    ("hpc-cap", "", cxxopts::value<int>())
    ("binary", "", cxxopts::value<bool>()->default_value("false"))
    ("version", "Print version information.")
    ("h,help", "Print this help.")
    ("l", "", cxxopts::value<float>())
    ("verbose", "", cxxopts::value<bool>()->default_value("false"));
  // clang-format on
}

void Configuration::parse(int argc, char **argv) {
  auto results = parser.parse(argc, argv);
  if (results.count("bam"))
    bam = results["bam"].as<std::string>();
  if (results.count("sfs"))
    sfs = results["sfs"].as<std::string>();
  if (results.count("poa"))
    poa = results["poa"].as<std::string>();
  if (results.count("clusters"))
    clusters = results["clusters"].as<std::string>();
  if (results.count("clip-clusters"))
    clip_clusters = results["clip-clusters"].as<std::string>();
  if (results.count("normal-contigs-bam"))
    normal_contigs_bam = results["normal-contigs-bam"].as<std::string>();
  if (results.count("index"))
    index = results["index"].as<std::string>();
  if (results.count("fastx"))
    fastq = results["fastx"]
                .as<std::string>(); // FIXME: use fastx here and in pingpong
  if (results.count("overlap"))
    overlap = results["overlap"].as<int>();
  if (results.count("bsize"))
    batch_size = results["bsize"].as<int>();
  if (results.count("omax"))
    max_output = results["omax"].as<int>();
  if (results.count("threads"))
    threads = results["threads"].as<int>();
  if (results.count("append"))
    append = results["append"].as<std::string>();
  if (results.count("reference"))
    reference = results["reference"].as<std::string>();
  if (results.count("min-sv-length"))
    // min_sv_length = max(25, results["min-sv-length"].as<int>());
    min_sv_length = results["min-sv-length"].as<int>();
  if (results.count("min-indel-length"))
    min_indel_length = results["min-indel-length"].as<int>();
  if (results.count("min-cluster-weight"))
    min_cluster_weight = results["min-cluster-weight"].as<int>();
  if (results.count("min-bnd-dist"))
    min_bnd_dist = results["min-bnd-dist"].as<int>();
  if (results.count("max-cluster-dist"))
    max_cluster_dist = results["max-cluster-dist"].as<int>();
  if (results.count("min-mapq"))
    min_mapq = results["min-mapq"].as<int>();
  if (results.count("accp"))
    accp = results["accp"].as<float>();
  if (results.count("germline-min-ro"))
    germline_min_ro = results["germline-min-ro"].as<float>();
  if (results.count("germline-min-reads"))
    germline_min_reads = results["germline-min-reads"].as<int>();
  if (results.count("germline-max-reads"))
    germline_max_reads = results["germline-max-reads"].as<int>();
  if (results.count("germline-diff-max"))
    germline_diff_max = results["germline-diff-max"].as<int>();
  if (results.count("germline-diff-min-len"))
    germline_diff_min_len = results["germline-diff-min-len"].as<int>();
  if (results.count("no-germline-diff") && results["no-germline-diff"].as<bool>())
    germline_diff = false;
  // diploid_poa, germline_realign and require_sfs_overlap now default to ON.
  // The old enabling flags stay accepted (no-ops) so existing command lines keep
  // working; the --no-* counterparts are the way to switch them off.
  if (results.count("no-diploid-poa") && results["no-diploid-poa"].as<bool>())
    diploid_poa = false;
  if (results.count("poa-min-freq"))
    poa_min_freq = results["poa-min-freq"].as<float>();
  if (results.count("no-germline-realign") &&
      results["no-germline-realign"].as<bool>())
    germline_realign = false;
  if (results.count("germline-realign-margin"))
    germline_realign_margin = results["germline-realign-margin"].as<int>();
  if (results.count("germline-realign-max-len"))
    germline_realign_max_len = results["germline-realign-max-len"].as<int>();
  if (results.count("l"))
    min_ratio = results["l"].as<float>();
  if (results.count("bed-exclusion")) {
    bed_exclusion = results["bed-exclusion"].as<std::string>();
    bed_filter.load(bed_exclusion);
  }
  binary = results["binary"].as<bool>();
  clipped = results["clipped"].as<bool>();
  clipped_fallback = !results["no-clipped-fallback"].as<bool>();
  if (results.count("no-require-sfs-overlap") &&
      results["no-require-sfs-overlap"].as<bool>())
    require_sfs_overlap = false;
  useht = !results["noht"].as<bool>();
  // noref = results["noref"].as<bool>();
  assemble = !(results["noassemble"].as<bool>());
  putative = !(results["noputative"].as<bool>());
  hpc = results["hpc"].as<bool>();
  if (results.count("hpc-cap"))
    hpc_cap = results["hpc-cap"].as<int>();
  version = results["version"].as<bool>();
  verbose = results["verbose"].as<bool>();
  help = results["help"].as<bool>();

  batch_size = (batch_size / threads) * threads;
}
