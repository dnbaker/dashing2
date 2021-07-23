#include "filterset.h"
namespace dashing2 {
    void FilterSet::finalize() {
        //reset(data_.begin(), data_.end(), bfexp_, k_ > 0 ? k_: -1);
        if(is_bf()) {
            const size_t nelem = data_.size();
            auto beg = data_.data(), end = &*data_.end();
            auto oldd = std::move(data_);
            const double dbits_per_el = std::log(1. / bfexp_) / std::pow(M_LN2, 2);
            int k = k_ > 0 ? k_: int(dbits_per_el * M_LN2);
            // Constant is std::pow(1. / std::log(2), 2.)
            size_t mem = nelem * dbits_per_el;
            mem = std::max((mem + mem % bits_per_reg) / bits_per_reg, size_t(1));
            k_ = k;
            std::fprintf(stderr, "nreg: %zu. k: %d. nbits: %zu\n", mem, k_, mem * 64);
            std::fprintf(stderr, "False positive rate %g using %d bits per el\n", bfexp_, k_);
            data_.resize(mem);
            std::fprintf(stderr, "data size: %zu. total mem: %zu\n", data_.size(), size_t(nelem * dbits_per_el));
            schism::Schismatic<uint64_t> mod_(mem * 64);
            // Maybe we batch this in the future?
            for(;beg != end; ++beg) {
                RNG mt(hash(*beg));
                for(int i = 0; i < k_; ++i) {
                    const uint64_t rem = mod_.mod(mt());
                    data_[(rem >> SHIFT)] |= (T(1) << (rem & MASK));
                }
            }
        } else {
            std::sort(data_.begin(), data_.end());
            //std::fprintf(stderr, "Sorted filterset of size %zu\n", data_.size());
            // Maybe replace with a faster sorter
        }
    }
} // namespace dashing2
