#include "bed_filter.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

void BedFilter::load(const string& bed_path) {
    if (bed_path.empty()) return;

    ifstream in(bed_path);
    if (!in.is_open()) {
        cerr << "Warning: Could not open BED file " << bed_path << endl;
        return;
    }

    map<string, vector<lib_interval_tree::interval<int>>> chr_intervals;

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        istringstream iss(line);
        string chrom;
        int s, e;
        if (iss >> chrom >> s >> e) {
            chr_intervals[chrom].push_back(lib_interval_tree::interval<int>(s, e - 1));
        }
    }

    for (auto& pair : chr_intervals) {
        lib_interval_tree::interval_tree_t<int> tree;
        for (const auto& ival : pair.second) {
            tree.insert(ival);
        }
        trees[pair.first] = move(tree);
    }
}

bool BedFilter::overlaps(const string& chrom, int s, int e) const {
    if (trees.empty()) return false;

    auto it = trees.find(chrom);
    if (it == trees.end()) return false;

    int low = min(s, e);
    int high = max(s, e);

    auto overlaps = it->second.overlap_find({low, high});
    return overlaps != it->second.end();
}
