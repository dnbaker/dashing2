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
#ifdef _OPENMP
        #pragma omp simd
#endif
        for(size_t i = 0; i < cs; ++i) {
            count_sketch_[i] += o.count_sketch_[i];
        }
    }
    return *this;
}

void Counter::reset() {
    std::fill(count_sketch_.begin(), count_sketch_.end(), 0.);
    if(!c64_.empty()) c64_.clear();
    if(!c128_.empty()) c128_.clear();
    if(!c64d_.empty()) c64d_.clear();
    if(!c128d_.empty()) c128d_.clear();
}
bool Counter::empty() const
{
    return c64_.empty() && c128_.empty() && std::find_if(count_sketch_.begin(), count_sketch_.end(), [](auto x) {return x != 0;}) == count_sketch_.end();
}

} // namespace dashing2
