#ifndef DASHING2_COUNTER_H__
#define DASHING2_COUNTER_H__
#include "enums.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include "sketch/div.h"
#include "hash.h"

namespace dashing2 {

struct Counter {
    ska::flat_hash_map<uint64_t, int32_t> c64_;
    ska::flat_hash_map<u128_t, int32_t, FHasher> c128_;
    ska::flat_hash_map<uint64_t, double> c64d_;
    ska::flat_hash_map<u128_t, double, FHasher> c128d_;
    std::vector<double> count_sketch_;
    schism::Schismatic<uint64_t> s64_;
    CountingType ct() const {return count_sketch_.size() ? COUNTSKETCH_COUNTING: EXACT_COUNTING;}
    Counter(size_t cssize=0): count_sketch_(cssize), s64_(cssize + !cssize) {}
    Counter &operator+=(const Counter &o);
    static constexpr uint64_t BM64 = 0x8000000000000000ull;
    void add(u128_t x) {
        if(ct() == EXACT_COUNTING) {
            auto it = c128_.find(x);
            if(it == c128_.end()) c128_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
            double inc = ct() == COUNTSKETCH_COUNTING && ((hv & BM64) == 0) ? -1.: 1.;
            count_sketch_[s64_.mod(hv)] += inc;
        }
    }
    void reset();
    bool empty() const {return c64_.empty() && c128_.empty() && std::find_if(count_sketch_.begin(), count_sketch_.end(), [](auto x) {return x != 0;}) == count_sketch_.end();}
    void add(u128_t x, double inc) {
        switch(ct()) {
            case COUNTSKETCH_COUNTING: case COUNTMIN_COUNTING: {
                const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
                const auto idx = s64_.mod(hv);
                if(ct() == COUNTSKETCH_COUNTING && ((hv & BM64) == 0)) inc = -inc;
                count_sketch_[idx] += inc;
            } break;
            case EXACT_COUNTING: {
                auto it = c128d_.find(x);
                if(it == c128d_.end())
                    c128d_.emplace(x, inc);
                else it->second += inc;
            } break;
            default: __builtin_unreachable();
        }
    }
    void add(uint64_t x, double inc) {
        switch(ct()) {
            case COUNTSKETCH_COUNTING: case COUNTMIN_COUNTING: {
                const auto hv = sketch::hash::WangHash::hash(x), idx = s64_.mod(hv);
                if(ct() == COUNTSKETCH_COUNTING && ((hv & BM64) == 0)) inc = -inc;
                count_sketch_[idx] += inc;
            } break;
            case EXACT_COUNTING: {
                auto it = c64d_.find(x);
                if(it == c64d_.end())
                    c64d_.emplace(x, inc);
                else it->second += inc;
            } break;
            default: __builtin_unreachable();
        }
    }
    void add(uint64_t x) {
        if(ct() == EXACT_COUNTING) {
            auto it = c64_.find(x);
            if(it == c64_.end()) c64_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(x);
            count_sketch_[s64_.mod(hv)] += (ct() == COUNTSKETCH_COUNTING && ((hv & BM64) == 0) ? -1.: 1.);
        }
    }
    template<typename IT, typename Alloc, typename CountT>
    void finalize(std::vector<IT, Alloc> &dst, std::vector<CountT> &countdst, const double threshold = 0) {
        std::vector<std::pair<IT, uint32_t>> tmp;
        if(ct() == EXACT_COUNTING) {
            auto update_if = [&](auto &src) {
                if(!src.empty()) {
                    tmp.resize(src.size());
                    auto tmpit = tmp.begin();
                    for(const auto &pair: src)
                        if(pair.second > threshold)
                            *tmpit++ = {pair.first, pair.second};
                    return true;
                }
                return false;
            };
            if(!update_if(c64_) && !update_if(c128_) && !update_if(c64d_) && !update_if(c128d_)) {
                std::fprintf(stderr, "Note: finalizing an empty vector\n");
                dst.clear(); countdst.clear();
                return;
            }
        } else {
            const size_t css = count_sketch_.size();
            for(size_t i = 0; i < css; ++i) {
                if(count_sketch_[i] > threshold)
                   tmp.push_back({i, count_sketch_[i]});
            }
        }
        std::sort(tmp.begin(), tmp.end()); // Sort the hashes
        dst.resize(tmp.size());
        countdst.resize(tmp.size());
        auto dstit = dst.begin();
        auto countdstit = countdst.begin();
        for(auto &pair: tmp) {
            *dstit++ = pair.first;
            *countdstit++ = pair.second;
        }
    }
    template<typename Sketch>
    void finalize(Sketch &dst, const double threshold = 0) {
        if(ct() == EXACT_COUNTING) {
            auto update_if = [&](auto &src) {
                if(!src.empty()) {
                    for(const auto &pair: src)
                        if(pair.second > threshold)
                            dst.update(pair.first, pair.second);
                    return true;
                }
                return false;
            };
            update_if(c64_) || update_if(c128_) || update_if(c64d_) || update_if(c128d_);
        } else {
            const size_t css = count_sketch_.size();
            OMP_ONLY(_Pragma("omp simd"))
            for(size_t i = 0; i < css; ++i)
                if(auto csv = std::abs(count_sketch_[i]); csv >= threshold)
                    dst.update(i, csv);
        }
    }
};

}

#endif

