#pragma once
#include "btree/map.h"
#include "enums.h"
#include "sketch/hash.h"
#include "sketch/div.h"
#include "sketch/common.h"

namespace dashing2 {

namespace hash = sketch::hash;
template<typename T, size_t pow2=false>
struct LazyOnePermSetSketch {
    size_t m_;
    // Solution: hash reversibly, track the maximum IDs
    
    std::vector<T> registers_;
    static_assert(std::is_integral_v<T> || std::is_same_v<T, u128_t>, "LazyOnePermSetSketch is to be used with integral types");
    std::vector<double> counts_;
    using SigT = std::conditional_t<(sizeof(T) == 4), float, std::conditional_t<(sizeof(T) == 8), double, long double>>;
    std::unique_ptr<std::vector<SigT>> as_sigs_;
    std::unique_ptr<std::vector<uint64_t>> original_ids_;
    std::unique_ptr<std::vector<uint32_t>> idcounts_;
    size_t mask_;
    int shift_;
    double count_threshold_;
    // MultiplyAddXoRot
    // is already enough to pass Rabbit/SmallCrush
    using Hasher =
#if 0
        hash::FusedReversible3<hash::XorMultiply, hash::RotL33, hash::MultiplyAddXoRot<31>>;
#else
        hash::MultiplyAddXoRot<31>;
#endif
    Hasher hasher_;
    schism::Schismatic<uint64_t> div_;
    double mincount_ = 0.;
    std::vector<btree::map<T, uint32_t>> potentials_;
    double card_ = -1.;
    LazyOnePermSetSketch(const LazyOnePermSetSketch &o) = default;
    LazyOnePermSetSketch& operator=(const LazyOnePermSetSketch &o) = default;
    LazyOnePermSetSketch& operator=(LazyOnePermSetSketch &&o) = default;
    LazyOnePermSetSketch(LazyOnePermSetSketch &&) = default;
    LazyOnePermSetSketch(size_t m, uint64_t seed=0x321b919a61cb41f7ul): hasher_(seed), div_(m) {
        if(pow2)
            m = sketch::integral::roundup(m);
        else if(m & 1) ++m;
        m_ = m;
        registers_.resize(m_, T(0));
        counts_.resize(m_);
        div_ = schism::Schismatic<uint64_t>(m_);
        mask_ = m_ - 1;
        shift_ = sketch::integral::ilog2(m_);
    }
    void set_mincount(double v) {
        if(v > 0.) {
            mincount_ = v;
            potentials_.resize(registers_.size());
        }
    }
    template<typename O>
    INLINE void update(T id,  O) {update(id);}
    size_t threshold(size_t idx) {
        auto mul = std::numeric_limits<T>::max() / 2 + 1;
        assert(mul * 2 == 0);
        //std::fprintf(stderr, "Expected %zu, found %zu. returning %zu\n", size_t(mul), 2 * mul, mul / (size() >> 1));
        return idx * (mul / (size() / 2));
    }
    INLINE void update(const T oid) {
        const T id = hasher_(oid);
        size_t idx;
        if constexpr(pow2) {
            idx = id & mask_;
            id >>= shift_;
        } else {
            auto di = div_.div(id);
            auto mo = id - m_ * di;
            assert(di == (id / m_));
            assert(mo == (id % m_));
            assert(di < id);
            idx = mo;
            //std::fprintf(stderr, "di = %zu, mo = %zu, idx = %zu, id = %zu\n", size_t(di), size_t(mo), size_t(idx), size_t(id));
        }
        assert(idx < size());
        auto &cref = counts_[idx];
        auto &rref = registers_[idx];
        if(mincount_ > 0.) {
            // If mincount > 0, then
            if(rref < id) {
                auto &pos = potentials_[idx];
                auto it = pos.find(id);
                if(it == pos.end()) it = pos.emplace(id, 1).first;
                else ++it->second;
                if(it->second >= mincount_) {
                    rref = id;
                    cref = it->second;
                    pos.erase(it);
                    for(auto pit = pos.begin(); pit != pos.end();) {
                        if(pit->first <= id) pit = pos.erase(pit);
                        else ++pit;
                    }
                    return;
                }
            } else if(rref == id) ++cref;
        } else if(rref < id) {
            rref = id; cref = 1.;
        } else if(rref == id) ++cref;
    }
    std::vector<uint32_t> &idcounts() {
        auto p = new std::vector<uint32_t>(size());
        idcounts_.reset(p);
        std::copy(counts_.begin(), counts_.end(), p->data());
        return *p;
    }
    auto &ids() {
        auto p = new std::vector<uint64_t>(registers_.size());
        original_ids_.reset(p);
        std::transform(registers_.begin(), registers_.end(), p->begin(), [this](T x) {return this->hasher_.inverse(x);});
        return *p;
    }
    SigT *data() {
        as_sigs_.reset(new std::vector<SigT>(registers_.size()));
        const auto modv = (std::numeric_limits<T>::max() / 2 + 1) / (size() / 2);
        size_t idx = 0;
        const SigT mul = -SigT(1) / size();
        const SigT omul = SigT(1) / modv;
        for(size_t i = 0; i < size(); ++i) {
            SigT v;
            const auto lv = registers_[i] / m_;
            if(lv) {
                v = mul * std::log(omul * lv);
            } else v = 0.;
            as_sigs_->operator[](i) = v;
        }
        return as_sigs_->data();
    }
    void reset() {
        std::memset(registers_.data(), 0, registers_.size() * sizeof(T));
        std::memset(counts_.data(), 0, counts_.size() * sizeof(double));
    }
    double getcard() {
        if(card_ > 0.) return card_;
        const auto modv = (std::numeric_limits<T>::max() / 2 + 1) / (size() / 2);
        long double inv = 1. / modv;
        long double sum = std::accumulate(registers_.begin(), registers_.end(), 0.L,
                    [modv,inv,m=m_](auto x, auto y) {
            return x + (modv - (y / m) % modv) * inv;
        });
        if(sum == 0.) return std::numeric_limits<double>::infinity();
        card_ = std::pow(size(), 2) / sum;
        return card_;
    }
    size_t size() const {return registers_.size();}
    template<typename T2=SigT>
    std::vector<T2> to_sigs() const {
        std::vector<T2> ret(size());
        std::transform(registers_.begin(), registers_.end(), ret.begin(), [sz2=size()/2](auto x) -> T2 {
            if(std::is_integral_v<T2>) {
                return x; // save as truncation/min hash value by twiddling
            } else {
                const auto modv = (std::numeric_limits<T>::max() / 2 + 1) / sz2;
                long double inv = 1. / modv;
                return -std::log((modv - (x / m) % modv) * inv);
            }
        });
        return ret;
    }
};

}
