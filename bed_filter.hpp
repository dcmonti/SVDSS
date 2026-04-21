#ifndef BED_FILTER_HPP
#define BED_FILTER_HPP

#include <string>
#include <vector>
#include <map>
#include "interval_tree.hpp"

using namespace std;

class BedFilter {
public:
    BedFilter() = default;
    
    void load(const string& bed_path);
    
    bool overlaps(const string& chrom, int s, int e) const;

private:
    map<string, lib_interval_tree::interval_tree_t<int>> trees;
};

#endif
