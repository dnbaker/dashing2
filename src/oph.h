#pragma once
#include "btree/map.h"
#include "enums.h"
#include "sketch/hash.h"
#include "sketch/div.h"
#include "sketch/common.h"

namespace dashing2 {

namespace hash = sketch::hash;
template<typename T, size_t pow2=false, typename Hasher = hash::MultiplyAddXoRot<31>>
struct LazyOnePermSetSketch {
private:
    static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>, "Must be integral and unsigned");
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
    uint64_t total_updates_ = 0;
    Hasher hasher_;
    // MultiplyAddXoRot
    // is already enough to pass Rabbit/SmallCrush
    schism::Schismatic<uint64_t> div_;
    double mincount_ = 0.;
    std::vector<btree::map<T, uint32_t>> potentials_;
    double card_ = -1.;
public:
    LazyOnePermSetSketch(const LazyOnePermSetSketch &o): div_(o.div_) {
        *this = o;
    }
    LazyOnePermSetSketch& operator=(const LazyOnePermSetSketch &o)
    {
        m_ = o.m_;
        registers_ = o.registers_;
        counts_ = o.counts_;
        mask_ = o.mask_;
        shift_ = o.shift_;
        count_threshold_ = o.count_threshold_;
        total_updates_ = o.total_updates_;
        div_ = o.div_;
        mincount_ = o.mincount_;
        potentials_ = o.potentials_;
        card_ = o.card_;
        if(o.as_sigs_) as_sigs_.reset(new std::decay_t<decltype(*as_sigs_)>(*o.as_sigs_));
        if(o.original_ids_) original_ids_.reset(new std::decay_t<decltype(*original_ids_)>(*o.original_ids_));
        if(o.idcounts_) idcounts_.reset(new std::decay_t<decltype(*idcounts_)>(*o.idcounts_));
        return *this;
    }
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
        reset();
    }
    void set_mincount(double v) {
        if(v > 1.) {
            mincount_ = v;
            potentials_.resize(size());
        }
    }
    template<typename O>
    INLINE void update(T id,  O) {update(id);}
    T decode(T x) const {
        return hasher_.inverse(x);
    }
    static INLINE constexpr T shasher_(T x) {
        return static_cast<T>(0xd63af43ad731df95) ^ ((x >> 31) ^ (x << 33));
    }
    INLINE void update(const T oid) {
        ++total_updates_;
        const T id = hasher_(oid);
        size_t idx = shasher_(id);
        if constexpr(pow2)
            idx = idx & mask_;
        else
            idx -= m_ * div_.div(idx);
        assert(idx < size());
        auto &cref = counts_[idx];
        auto &rref = registers_[idx];
        if(mincount_ > 1.) {
            // If mincount > 0, then
            if(rref > id) {
                auto &pos = potentials_[idx];
                auto it = pos.find(id);
                if(it == pos.end()) it = pos.emplace(id, 1).first;
                else ++it->second;
                if(it->second >= mincount_) {
                    rref = id;
                    cref = it->second;
                    for(auto pit = pos.begin(); pit != pos.end();) {
                        if(pit->first >= id)
                        pit = pos.erase(pit);
                        else ++pit;
                    }
                    return;
                }
            } else cref += (rref == id);
        } else {
            if(rref > id) {
                rref = id; cref = 1.;
            } else cref += (rref == id);
        }
    }

    static constexpr long double omul =
        sizeof(T) == 16 ? 0x1p-128L:
        sizeof(T) == 8 ? 0x1p-64L:
        sizeof(T) == 4 ? 0x1p-32L:
        sizeof(T) == 2 ? 0x1p-16L:
        sizeof(T) == 1 ? 0x1p-8L: 0.L;
    static_assert(omul != 0.L, "sanity check");
    template<typename T2=SigT>
    std::vector<T2> to_sigs() const {
        std::vector<T2> ret(size());
        std::transform(registers_.begin(), registers_.end(), ret.begin(), [sz2=size()/2](auto x) -> T2 {
            if(std::is_integral_v<T2>) {
                return x; // save as truncation/min hash value by twiddling
            } else {
                return -std::log((get_modv() - x) * omul);
            }
        });
        return ret;
    }
    void reset() {
        std::fill_n(registers_.data(), registers_.size(), std::numeric_limits<T>::max());
        std::memset(counts_.data(), 0, counts_.size() * sizeof(double));
        as_sigs_.reset();
        card_ = -1.;
        for(auto &p: potentials_) p.clear();
    }
    double getcard() {
        if(card_ > 0.) return card_;
        long double sum = std::accumulate(registers_.begin(), registers_.end(), 0.L,
            [](auto x, auto y) {return x + y * omul;}
        );
        if(!sum) return std::numeric_limits<double>::infinity();
        return m_ * (m_ / sum);
    }
    SigT *data() {
        if(as_sigs_) return as_sigs_->data();
        as_sigs_.reset(new std::vector<SigT>(registers_.size()));
        auto asp = as_sigs_->data();
        const long double mul = -SigT(1) / (m_ - std::count(registers_.begin(), registers_.end(), std::numeric_limits<T>::max()));
        for(size_t i = 0; i < m_; ++i) {
            if(registers_[i] == std::numeric_limits<T>::max()) continue;
            asp[i] = mul * std::log(omul * (std::numeric_limits<T>::max() - registers_[i] + 1));
        }
        return asp;
    }
    std::vector<uint64_t> &ids() {
        auto p = new std::vector<uint64_t>(registers_.size());
        original_ids_.reset(p);
        std::transform(registers_.begin(), registers_.end(), p->begin(), [this](T x) {return this->decode(x);});
        return *p;
    }
    std::vector<uint32_t> &idcounts() {
        auto p = new std::vector<uint32_t>(size());
        idcounts_.reset(p);
        std::copy(counts_.begin(), counts_.end(), p->data());
        return *p;
    }
    static constexpr size_t get_modv() {return std::numeric_limits<T>::max();}
    size_t total_updates() const {return total_updates_;}
    size_t size() const {return m_;}
};
} // namespace dashing2
