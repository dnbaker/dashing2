#ifndef DASHING2_INDEX_BUILD_H__
#define DASHING2_INDEX_BUILD_H__
#include "sketch/ssi.h"
#include "src/cmp_main.h"
namespace dashing2 {

using PairT = std::pair<LSHDistType, LSHIDType>;

std::vector<std::vector<PairT>>
build_index(SetSketchIndex<uint64_t, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result,
    const bool index_compressed=false);

} // namespace dashing2

#endif
