#include "counter.h"
#include "merge.h"
namespace dashing2 {

Counter &Counter::operator+=(const Counter &o) {
    if(ct() != o.ct() ||count_sketch_.size() != o.count_sketch_.size())
        throw std::invalid_argument("Counters do not share parameters");
    if(ct() == EXACT_COUNTING) {
        if(c64_.size()) {
            merge(c64_, o.c64_);
        } else if(c128_.size()) {
            merge(c128_, o.c128_);
        } else if(c64d_.size()) {
            merge(c64d_, o.c64d_);
        } else if(c128d_.size()) {
            merge(c128d_, o.c128d_);
        }
    } else {
        const size_t cs = count_sketch_.size();
        #pragma omp simd
        for(size_t i = 0; i < cs; ++i) {
            count_sketch_[i] += o.count_sketch_[i];
        }
    }
    return *this;
}

}
