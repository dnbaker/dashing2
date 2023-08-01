#ifndef SKETCH_CORE_H__
#define SKETCH_CORE_H__
#include "fastxsketch.h"
#include "lfsketch.h"
#include "bwsketch.h"
#include "bedsketch.h"
namespace dashing2 {
std::vector<std::pair<size_t, uint64_t>> get_filesizes(const std::vector<std::string> &paths);
SketchingResult &sketch_core(SketchingResult &, Dashing2DistOptions &opts, const std::vector<std::string> &paths, std::string &outfile);
}

#endif
