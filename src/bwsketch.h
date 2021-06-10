#ifndef DASHING2_BWSKTCH_H__
#define DASHING2_BWSKTCH_H__
#include "d2.h"
namespace dashing2 {

using BigWigSketchResult = IntervalSketchResult;

BigWigSketchResult bw2sketch(std::string path, const Dashing2Options &opts);

}
#endif
