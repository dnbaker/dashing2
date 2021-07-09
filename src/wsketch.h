#ifndef SIMPLE_MH_H__
#define SIMPLE_MH_H__
#include <utility>
#include <cstdint>
#include <vector>
#include "d2.h"
namespace dashing2 {
using std::uint64_t;


struct SimpleMHRet: public std::tuple<std::vector<RegT>, std::vector<uint64_t>, std::vector<uint64_t>, double> {
    using Tup =  std::tuple<std::vector<RegT>, std::vector<uint64_t>, std::vector<uint64_t>, double>;
    Tup &tup() {return *static_cast<Tup *>(this);}
    auto &sigs() {return std::get<0>(*this);}
    auto &hashes() {return std::get<1>(*this);}
    auto &ids() {return std::get<2>(*this);}
    auto &total_weight() {return std::get<3>(*this);}
    SimpleMHRet &operator=(Tup &&o) {
        Tup::operator=(o);
        return *this;
    }
    SimpleMHRet &operator=(const Tup &o) {
        Tup::operator=(o);
        return *this;
    }
};

#if 0
SimpleMHRet minhashf64u64(const double *, const uint64_t *, size_t, size_t, int usepmh);
SimpleMHRet minhashf64u32(const double *, const uint32_t *, size_t, size_t, int usepmh);
SimpleMHRet minhashf32u64(const float *, const uint64_t *, size_t, size_t, int usepmh);
SimpleMHRet minhashf32u32(const float *, const uint32_t *, size_t, size_t, int usepmh);
static inline SimpleMHRet minhash(const double *w, const uint64_t *ids, size_t n, size_t m, int usepmh) {return minhashf64u64(w, ids, n, m, usepmh);}
static inline SimpleMHRet minhash(const double *w, const uint32_t *ids, size_t n, size_t m, int usepmh) {return minhashf64u32(w, ids, n, m, usepmh);}
static inline SimpleMHRet minhash(const float *w, const uint64_t *ids, size_t n, size_t m, int usepmh) {return minhashf32u64(w, ids, n, m, usepmh);}
static inline SimpleMHRet minhash(const float *w, const uint32_t *ids, size_t n, size_t m, int usepmh) {return minhashf32u32(w, ids, n, m, usepmh);}
#endif

SimpleMHRet wmh_from_file(std::string idpath, std::string cpath, size_t sksz, bool usepmh, bool usef32=false, bool wordids=false);
}


#endif
