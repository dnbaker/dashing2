#ifndef DASHING2_FILTERSET_H__
#define DASHING2_FILTERSET_H__
#include "sketch/sseutil.h"
#include "sketch/div.h"
#include "sketch/macros.h"
#include "aesctr/wy.h"
#include "aligned_vector.h"
#include "interpbound.h"
#include <climits>
#include <cmath>


namespace dashing2 {
static constexpr int ceilog2(size_t n) {
    if(n == 0 || (n & (n - 1)) == 0) {
        switch(n) {
            case 1: return 0;
            case 2: return 1;
            case 4: return 2;
            case 8: return 3;
            case 16: return 4;
            case 32: return 5;
            case 64: return 6;
            case 128: return 7;
            case 256: return 8;
        }
    }
    return -1;
}


class FilterSet {
    using T = uint64_t;
    aligned::vector<T> data_;
    double bfexp_;
    int k_ = 0;
    //  Table size is
#ifndef M_LN2
    static constexpr double M_LN2 = 0.6931471805599453;
#endif
    using RNG = wy::WyRand<uint64_t, 2>;
    static constexpr size_t bits_per_reg = sizeof(T) * CHAR_BIT;
public:
    std::string to_string() const {
        if(is_bf()) {
            return std::string("FilterSetBloomFilter-size=") + std::to_string(data_.size()) + ",k=" + std::to_string(k_) + ",err=" + std::to_string(bfexp_);
        }
        return std::string("FilterSetSortedHashSet-size=") + std::to_string(data_.size());
    }
    template<typename IT>
    void reset(IT beg, IT end, double bfexp=-1., int ktouse=-1) {
        *this = FilterSet(beg, end, bfexp, ktouse);
    }
    FilterSet &operator=(const FilterSet &o) = default;
    FilterSet &operator=(FilterSet &&o){
        data_ = std::move(o.data_);
        bfexp_ = o.bfexp_;
        k_ = o.k_;
        return *this;
    }
    template<typename IT>
    FilterSet(IT beg, IT end, double bfexp=-1., int ktouse=-1): bfexp_(bfexp) {
        std::fprintf(stderr, "bfexp: %g\n", bfexp_);
        const size_t nelem = std::distance(beg, end);
        if(is_bf()) {
            const double dbits_per_el = std::log(1. / bfexp) / std::pow(M_LN2, 2);
            int k = ktouse > 0 ? ktouse: int(dbits_per_el * M_LN2);
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
                RNG mt(*beg);
                for(int i = 0; i < k_; ++i) {
                    auto v = mt();
                    const size_t rem = mod_.mod(v);
                    auto &reg = data_[(rem >> SHIFT)];
                    reg |= (T(1) << (rem & MASK));
                    data_[(rem >> SHIFT)] |= (T(1) << (rem & MASK));
                }
            }
        } else {
            data_.reserve(nelem);
            std::copy(beg, end, std::back_inserter(data_));
            std::sort(data_.begin(), data_.end());
            // Maybe replace with a faster sorter
        }
    }
    bool is_bf() const {return bfexp_ > 0.;}
    size_t tablesize() const {return (data_.size());}
    static constexpr int SHIFT = ceilog2(sizeof(T) * CHAR_BIT);
    static constexpr uint64_t MASK = (size_t(1) << SHIFT) - 1;
    size_t nbits() const {return tablesize() << SHIFT;}
    bool in_set(T x) const {
        if(is_bf()) {
            // TODO: batch to sort for efficiency
            RNG rng(x);
            schism::Schismatic<uint64_t> mod_(nbits());
            // Maybe this can even be SIMDIfied?
            for(int i = 0; i < k_; ++i) {
                const size_t v = rng();
                const size_t rem = mod_.mod(v);
                //std::fprintf(stderr, "Getting hash %zu from rem %zu - %d at index %zu at bit %d\n", v, rem, i, (rem >> SHIFT), (rem & MASK));
                if((data_[rem >> SHIFT] & (uint64_t(1) << (rem & MASK))) == 0) {
                    //std::fprintf(stderr, "Data %zu at bit %d is 0\n", rem >> SHIFT, (rem & MASK));
                    return false;
                }
            }
            return true;
        }
        //std::fprintf(stderr, "Testing with hash set lookup, and %s\n", std::is_sorted(data_.begin(), data_.end()) ? "is sorted": "is not sorted");
        // Replace with interpolation search in the future
        const T *p = data_.data(),
                *e = p + data_.size(),
                *it = std::lower_bound(p, e, x);
        // Equivalent, but no logging
        return it != e && *it == x;
    }
};

}

#endif /* DASHING2_FILTERSET_H__ */
