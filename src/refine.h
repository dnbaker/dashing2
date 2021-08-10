#pragma once
#ifndef DASHING2_REFINE_H__
#define DASHING2_REFINE_H__
#include "index_build.h"

namespace dashing2 {
void refine_results(std::vector<pqueue> &lists, const Dashing2DistOptions &opts, const SketchingResult &result);
}

#endif
