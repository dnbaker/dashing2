#pragma once
#include "btree/map.h"
#include "enums.h"
#include "sketch/div.h"
#include "sketch/common.h"
#include "hash.h"

namespace dashing2 {

//using u128_t = __uint128_t;

template<uint64_t base_seed>
struct CEIAdd {
    static constexpr uint64_t seed_ = base_seed;
    static constexpr uint32_t seed32_ = uint32_t(base_seed);
    static constexpr uint64_t seed2_ = sketch::hash::WangHash::hash(base_seed);
    static constexpr __uint128_t sl = (__uint128_t(seed_) << 64) | seed2_;
    template<typename...Args>
    CEIAdd(Args &&...) {}

    INLINE uint64_t inverse(uint64_t hv) const {
        return hv - seed_;
    }
    INLINE uint32_t inverse(uint32_t hv) const { return hv - seed32_;}
    INLINE uint64_t operator()(uint64_t h) const {return h + seed_;}
    INLINE uint32_t operator()(uint32_t h) const {return h + seed32_;}
    INLINE __uint128_t inverse(__uint128_t hv) const {
        return hv - sl;
    }
    INLINE __uint128_t operator()(__uint128_t hv) const {
        return hv + sl;
    }
};
using sketch::hash::CEIXOR;
using sketch::hash::CEIMul;
#if 0
// USE_SIMPLE_REVHASH
using BHasher = sketch::hash::CEIFused<CEIXOR<0x533f8c2151b20f97>, CEIMul<0x9a98567ed20c127d>, CEIXOR<0x691a9d706391077a>>;
#else
using BHasher = sketch::hash::WangHash;
#endif
struct DHasher {
    BHasher bh;
    uint64_t seed_;
    u128_t seed2_;
    DHasher(uint64_t x): seed_(std::mt19937_64(x)()), seed2_((u128_t(sketch::hash::WangHash::hash(x)) << 64) | seed_)
    {
#ifndef NDEBUG
        const uint64_t tmpv = 133348;
        assert(bh.inverse(bh(tmpv)) == tmpv);
        assert(inverse(operator()(tmpv)) == tmpv);
#endif
    }
    //uint32_t operator()(uint32_t x) const {return bh(x ^ static_cast<uint32_t>(seed_));}
    //uint32_t operator()(int32_t x) const {return operator()(static_cast<uint32_t>(x));}
    template<typename IT, typename=std::enable_if_t<std::is_integral_v<IT> && !std::is_same_v<IT, uint64_t> && (sizeof(IT) <= 8)>>
    uint64_t operator()(IT x) const {return bh(uint64_t(x) ^ seed_);}
    uint64_t operator()(uint64_t x) const {return bh(x ^ seed_);}
    uint64_t operator()(int32_t x) const {return bh(x ^ seed_);}
    u128_t operator()(u128_t x) const {return bh(x ^ seed2_);}
#if 0
    uint32_t inverse(uint32_t x) const {
        return static_cast<uint32_t>(seed_) ^ bh.inverse(x);
    }
    uint32_t inverse(int32_t x) const {return inverse(static_cast<uint32_t>(x));}
    uint64_t inverse(int64_t x) const {return inverse(static_cast<uint64_t>(x));}
#endif
    uint64_t inverse(uint64_t x) const {
        return seed_ ^ bh.inverse(x);
    }
    u128_t inverse(u128_t x) const {
        return seed2_ ^ bh.inverse(x);
    }
};



namespace hash = sketch::hash;
template<typename T, size_t pow2=false, typename Hasher = DHasher>
struct LazyOnePermSetSketch {
private:
    static_assert((std::is_integral_v<T> && std::is_unsigned_v<T>) || std::is_same_v<T, u128_t>, "Must be integral and unsigned");
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
    std::vector<flat_hash_map<T, uint32_t>> potentials_;
    double card_ = -1.;
public:
    LazyOnePermSetSketch(const LazyOnePermSetSketch &o): hasher_(o.hasher_), div_(o.div_) {
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
        registers_.resize(m_, T(-1));
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
    static INLINE constexpr uint32_t shasher_(uint32_t x) {
        return static_cast<uint32_t>(0xd63af43a) ^ ((x >> 15) ^ (x << 17));
    }
    static INLINE constexpr uint64_t shasher_(uint64_t x) {
        return static_cast<T>(0xd63af43ad731df95) ^ ((x >> 31) ^ (x << 33));
    }
    template<typename OT> static INLINE constexpr OT shasher_(OT x) {
        if(sizeof(OT) == 4) {return shasher_(uint32_t(x));}
        else if(sizeof(OT) == 16) {return shasher_(uint64_t(x) ^ shasher_(uint64_t(x >> 64)));}
        return shasher_(uint64_t(x));
    }
    INLINE void update(const T oid) {
        ++total_updates_;
        const T id = hasher_(oid);
        size_t idx = shasher_(id);
        if constexpr(pow2)
            idx = idx & mask_;
        else
            idx -= m_ * div_.div(idx);
        std::fprintf(stderr, "Pos: %zu. . Found pos: %zu. hash value: %zu. finalized-hash value %zu.\n", shasher_(id) % mask_, idx, id, shasher_(id));
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
        std::fill_n(registers_.data(), registers_.size(), T(-1));
        std::memset(counts_.data(), 0, counts_.size() * sizeof(double));
        as_sigs_.reset();
        card_ = -1.;
        for(auto &p: potentials_) p.clear();
        total_updates_ = 0;
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
        std::transform(registers_.begin(), registers_.end(), p->begin(), [this](const T &x) -> T {
            return this->decode(x);
        });
        return *p;
    }
    std::vector<uint32_t> &idcounts() {
        auto p = new std::vector<uint32_t>(size());
        idcounts_.reset(p);
        std::copy(counts_.begin(), counts_.end(), p->data());
        return *p;
    }
    static constexpr T get_modv() {return std::numeric_limits<T>::max();}
    size_t total_updates() const {return total_updates_;}
    size_t size() const {return m_;}
};
} // namespace dashing2
