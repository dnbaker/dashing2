#ifndef SIMPLE_MH_H__
#define SIMPLE_MH_H__
#include <utility>
#include <cstdint>
#include <vector>
#include "d2.h"
namespace dashing2 {
using std::uint64_t;


using SimpleMHRet = std::tuple<std::vector<RegT>, std::vector<uint64_t>, std::vector<uint64_t>, double>;
#pragma message("SimpleMHRet")

SimpleMHRet minhashf64u64(const double *, const uint64_t *, size_t, size_t, bool usepmh);
SimpleMHRet minhashf64u32(const double *, const uint32_t *, size_t, size_t, bool usepmh);
SimpleMHRet minhashf32u64(const float *, const uint64_t *, size_t, size_t, bool usepmh);
SimpleMHRet minhashf32u32(const float *, const uint32_t *, size_t, size_t, bool usepmh);
SimpleMHRet minhash(const double *w, const uint64_t *ids, size_t n, size_t m, bool usepmh) {return minhashf64u64(w, ids, n, m, usepmh);}
SimpleMHRet minhash(const double *w, const uint32_t *ids, size_t n, size_t m, bool usepmh) {return minhashf64u32(w, ids, n, m, usepmh);}
SimpleMHRet minhash(const float *w, const uint64_t *ids, size_t n, size_t m, bool usepmh) {return minhashf32u64(w, ids, n, m, usepmh);}
SimpleMHRet minhash(const float *w, const uint32_t *ids, size_t n, size_t m, bool usepmh) {return minhashf32u32(w, ids, n, m, usepmh);}

std::pair<std::vector<RegT>, std::vector<uint64_t>> wmh_from_file(std::string path, size_t sksz, bool usepmh);
}


#endif
