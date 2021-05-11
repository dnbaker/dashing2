#ifndef DASHING2_COUNTER_H__
#define DASHING2_COUNTER_H__
#include "enums.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include "sketch/div.h"
#include "hash.h"

namespace dashing2 {

struct Counter {
    CountingType ct_;
    ska::flat_hash_map<uint64_t, uint32_t> c64_;
    ska::flat_hash_map<u128_t, uint32_t, FHasher> c128_;
    ska::flat_hash_map<uint64_t, double> c64d_;
    ska::flat_hash_map<u128_t, double, FHasher> c128d_;
    std::vector<double> count_sketch_;
    schism::Schismatic<uint64_t> s64_;
    Counter(size_t cssize=0): ct_(cssize ? COUNTSKETCH_COUNTING: EXACT_COUNTING), count_sketch_(cssize), s64_(cssize + !cssize) {}
    Counter &operator+=(const Counter &o);
    static constexpr uint64_t BM64 = 0x8000000000000000ull;
    void add(u128_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c128_.find(x);
            if(it == c128_.end()) c128_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
    bool empty() const {return c64_.empty() && c128_.empty() && std::find_if(count_sketch_.begin(), count_sketch_.end(), [](auto x) {return x != 0;}) == count_sketch_.end();}
    void add(u128_t x, double inc) {
        if(ct_ == COUNTSKETCH_COUNTING) {
            const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1.: -1.) * inc;
        } else {
            auto it = c128d_.find(x);
            if(it == c128d_.end())
                c128d_.emplace(x, inc);
            else it->second += inc;
        }
    }
    void add(uint64_t x, double inc) {
        if(ct_ == COUNTSKETCH_COUNTING) {
            const auto hv = sketch::hash::WangHash::hash(x);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1.: -1.) * inc;
        } else {
            auto it = c64d_.find(x);
            if(it == c64d_.end())
                c64d_.emplace(x, inc);
            else it->second += inc;
        }
    }
    void add(uint64_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c64_.find(x);
            if(it == c64_.end()) c64_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(x);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
};

}

#endif

