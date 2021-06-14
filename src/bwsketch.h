#ifndef DASHING2_BWSKTCH_H__
#define DASHING2_BWSKTCH_H__
#include "d2.h"
namespace dashing2 {

using BigWigSketchResult = IntervalSketchResult;

BigWigSketchResult bw2sketch(std::string path, const Dashing2Options &opts);
std::vector<RegT> reduce(const ska::flat_hash_map<std::string, std::vector<RegT>> &map);
#ifndef BLOCKS_PER_ITER
#define BLOCKS_PER_ITER 4000000
#endif
static constexpr size_t blocks_per_iter = BLOCKS_PER_ITER;

}
#endif
