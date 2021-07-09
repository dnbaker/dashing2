#ifndef DASHING2_FILTERSET_H__
#define DASHING2_FILTERSET_H__
#include "sketch/sseutil.h"
#include "sketch/hash.h"
#include "sketch/div.h"
#include "sketch/macros.h"
#include "aesctr/wy.h"
#include "aligned_vector.h"
#include "interpbound.h"
#include "enums.h"
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
using namespace sketch::hash;


class FilterSet {
    using T = uint64_t;
    aligned::vector<T> data_;
    double bfexp_;
    int k_ = 0;
    //  Table size is
#ifndef M_LN2
    static constexpr double M_LN2 = 0.6931471805599453;
#endif
    using RNG = wy::WyRand<uint64_t, 0>;
    CEIFused<CEIXOR<0x533f8c2151b20f97>, CEIMul<0x9a98567ed20c127d>, CEIXOR<0x691a9d706391077a>>
        fshasher_;
    static constexpr size_t bits_per_reg = sizeof(T) * CHAR_BIT;
    uint64_t hash(uint64_t x) const {return fshasher_(x);}
    uint64_t hash(u128_t x) const {return fshasher_(uint64_t(x)) - fshasher_(uint64_t(x>>64));}
    template<typename T>
    uint64_t hash(T x) {
        if constexpr(sizeof(T) <= 8) return hash(uint64_t(x));
        else return hash(u128_t(x));
    }
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
    FilterSet(double bfexp=-1., int k=-1): bfexp_(bfexp), k_(k) {
    }
    void add(u128_t item) {
        data_.push_back(item >> 64);
        data_.push_back(item);
    }
    void add(uint64_t item) {
        data_.push_back(item);
    }
    void finalize() {
        reset(data_.begin(), data_.end(), bfexp_, k_ > 0 ? k_: -1);
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
                RNG mt(hash(*beg));
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
            if(&*beg != &data_[0])
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
#if 0
private:
    std::vector<std::tuple<uint64_t, uint64_t, RNG>> tmpvs;
    std::mutex mut;
public:
    std::vector<uint64_t> in_set(uint64_t *beg, uint64_t *end) {
        std::lock_guard<std::mutex> lock(mut);
        std::ptrdiff_t n = end - beg;
        if(n == 0) return {};
        if(n < 0) throw std::invalid_argument("Iterators are not ordered correctly.");
        std::vector<uint64_t> ret(((end - beg) + 63) / 64);
        if(is_bf()) {
            if(tmpvs.capacity() < size_t(n)) tmpvs.resize(n);
            for(ptrdiff_t i = 0; i < n; ++i) {
                const auto v = beg[i];
                std::get<1>(tmpvs[i]) = i;
                std::get<2>(tmpvs[i]) = RNG(v);
                std::get<0>(tmpvs[i]) = i;
            }
            std::sort(tmpvs.begin(), tmpvs.end());
            auto tmpp = tmpids.data();
            auto tvpp = tmpvs.data();
            std::transform(beg, end, rngs.begin(), RNG);
            std::transform(rngs.begin(), rngs.end(), tmpvs.begin(), [mp=&mod_](auto x) {return mp->mod(x);});
            std::transform(beg, end, tmpvs.data(), [mp=&mod_](auto x) {return mp->mod(RNG(x)());});
            std::iota(tmpp, tmpp + n, size_t(0));
            std::sort(tmpp, tmpp, [&](auto x, auto y) {return tmpvs[x] < tmpvs[y];});
            for(size_t i = 0; i < n; ++i) {
            }
            std::sort(&tmpvs[0], &tmpvs[n]);
        } else {
        }
        return ret;
    }
#endif
    bool in_set(u128_t x) const {
        if(is_bf()) {
            // TODO: batch to sort for efficiency -- see above drafting
            RNG rng((x >> 64) - (x & u128_t(0xFFFFFFFFFFFFFFFF)));
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
        const  u128_t *p = (const u128_t *)data_.data(), *e = p + (data_.size() / 2), *it = std::lower_bound(p, e, x);
        return it != e && *it == x;
    }
    bool in_set(T x) const {
        if(is_bf()) {
            // TODO: batch to sort for efficiency -- see above drafting
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
        const T *p = data_.data(), *e = p + data_.size(), *it = std::lower_bound(p, e, x);
        // Equivalent, but no logging
        return it != e && *it == x;
    }
};

FilterSet from_fastx(const std::string &path);

}

#endif /* DASHING2_FILTERSET_H__ */
