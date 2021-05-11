#ifndef DASHING2_BEDSKTCH_H__
#define DASHING2_BEDSKTCH_H__
#include "d2.h"

ska::flat_hash_map<std::string, std::vector<RegT>> bw2sketch(std::string path, const ParseOptions &opts);
std::vector<RegT> reduce(const ska::flat_hash_map<std::string, std::vector<RegT>> &map);
#endif
