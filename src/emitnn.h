#pragma once
#ifndef DASHING2_EMIT_NN_H__
#define DASHING2_EMIT_NN_H__
#include "index_build.h"

namespace dashing2 {
void emit_neighbors(std::vector<pqueue> &lists, const Dashing2DistOptions &opts, const SketchingResult &result);
}

#endif
