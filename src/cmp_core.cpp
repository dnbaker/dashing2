#include "cmp_main.h"
#include "sketch/count_eq.h"
#include "sketch/hash.h"
#include "index_build.h"
#include "refine.h"
#include "emitnn.h"
#include "mio.hpp"
#include "wcompare.h"
#include "levenshtein-sse.hpp"
#include "options.h"


namespace dashing2 {
std::pair<std::vector<LSHIDType>, std::vector<std::vector<LSHIDType>>> dedup_core(sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, const Dashing2DistOptions &opts, const SketchingResult &result);
void dedup_emit(const std::vector<LSHIDType> &, const std::vector<std::vector<LSHIDType>> &constituents, const Dashing2DistOptions &opts, const SketchingResult &result);

static INLINE uint64_t reg2sig(__float128 x) {
    const uint64_t *const p = (const uint64_t *)&x;
    return sketch::hash::WangHash::hash(sketch::hash::WangHash::hash(p[0] ^ 0xa3407fb23cd20eful) ^ p[1]);
}
static INLINE uint64_t reg2sig(long double x) {
    const uint64_t *const p = (const uint64_t *)&x;
    return sketch::hash::WangHash::hash(sketch::hash::WangHash::hash(p[0] ^ 0xa3407fb23cd20eful) ^ p[1]);
}
static INLINE uint64_t reg2sig(double x) {
    uint64_t v;
    std::memcpy(&v, &x, 8);
    return sketch::hash::WangHash::hash(v ^ 0xa3407fb23cd20eful);
}
static INLINE uint64_t reg2sig(float x) {
    uint32_t v;
    std::memcpy(&v, &x, 4);
    return sketch::hash::WangHash::hash(v ^ 0xa3407fb23cd20eful);
}

template<typename T>
size_t simdcount(const T *ptr, size_t n, const T v=static_cast<T>(0)) {
    size_t ret = 0;
    #pragma omp simd reduction(+:ret)
    for(size_t i = 0; i < n; ++i) {
        ret += ptr[i] == v;
    }
    return ret;
}

namespace detail{
#if 0
INLINE void setnthbit(uint64_t *ptr, size_t index, bool val) {
    ptr[index / 64] |= uint64_t(val) << (index % 64);
}
template<typename T> INLINE void setnthbit(T *ptr, size_t index, bool val) {
    return setnthbit(reinterpret_cast<uint64_t *>(ptr), index, val);
}

static INLINE uint64_t getnthbit(const uint64_t *ptr, size_t index) {
    return (ptr[index / 64] >> (index % 64)) & 1u;
}

template<typename T> INLINE T getnthbit(const T *ptr, size_t index) {
    return T(getnthbit(reinterpret_cast<const uint8_t *>(ptr), index));
}
INLINE uint64_t getnthbit(uint64_t val, size_t index) {
    return getnthbit(&val, index);
}


#if __SSE2__

INLINE bool is_all_zeros(__m128i x) {
#if __SSE4_1__
    return _mm_test_all_zeros(x, x);
#else
    return _mm_movemask_epi8(_mm_cmpeq_epi32(x,_mm_setzero_si128())) == 0xFFFF;
#endif
}

INLINE uint64_t matching_bits(const __m128i *s1, const __m128i *s2, uint16_t b) {
    --b;
    __m128i match = ~(*s1++ ^ *s2++);
    while(b >= 4) {
        match &= (~(*s1 ^ *s2)) & (~(s1[1] ^ s2[1])) & (~(s1[2] ^ s2[2])) & (~(s1[3] ^ s2[3]));
        s1 += 4;
        s2 += 4;
        b -= 4;
        if(is_all_zeros(match)) return 0;
    }
    while(b-- && !is_all_zeros(match)) match &= ~(*s1++ ^ *s2++);
    return popcount(common::vatpos(match, 0)) + popcount(common::vatpos(match, 1));
}
#endif

#if __AVX2__
INLINE bool is_all_zeros(__m256i x) {return _mm256_testz_si256(x, x);}
INLINE auto matching_bits(const __m256i *s1, const __m256i *s2, uint16_t b) {
    // Only do popcnt if match is nonzero
    --b;
    __m256i match = ~(*s1++ ^ *s2++);
    while(b >= 4) {
        match &= (~(*s1 ^ *s2)) & (~(s1[1] ^ s2[1])) & (~(s1[2] ^ s2[2])) & (~(s1[3] ^ s2[3]));
        b -= 4;
        s1 += 4;
        s2 += 4;
        if(is_all_zeros(match)) return match;
    }
    switch(b) {
        case 3: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 2: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 1: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 0: default: ;
    }
    return is_all_zeros(match) ? match: popcnt256(match);
}
#endif


#if HAS_AVX_512
INLINE bool is_all_zeros(__m512i x) {
    return _mm512_test_epi64_mask(x, x) == 0;
}
INLINE auto sbit_accum(const __m512i *s1, const __m512i *s2, uint16_t b) {
    __m512i match = ~(*s1++ ^ *s2++);
    while(--b)
        match &= ~(*s1++ ^ *s2++);
    return match;
}
INLINE auto matching_bits(const __m512i *s1, const __m512i *s2, uint16_t b) {
#define ONE_ITER do {\
        match &= (~(*s1 ^ *s2)) & (~(s1[1] ^ s2[1])) & (~(s1[2] ^ s2[2])) & (~(s1[3] ^ s2[3])); \
        if(is_all_zeros(match)) return match;\
        b -= 4; s1 += 4; s2 += 4;\
        } while(0)
    // Only do popcnt if match is nonzero
    --b;
    __m512i match = ~(*s1++ ^ *s2++);
    while(b >= 4) {
        switch(b >> 2) {
            case 8: ONE_ITER;[[fallthrough]];
            case 7: ONE_ITER;[[fallthrough]];
            case 6: ONE_ITER;[[fallthrough]];
            case 5: ONE_ITER;[[fallthrough]];
            case 4: ONE_ITER;[[fallthrough]];
            case 3: ONE_ITER;[[fallthrough]];
            case 2: ONE_ITER;[[fallthrough]];
            case 1: ONE_ITER;
        }
    }
#undef ONE_ITER
    switch(b) {
        case 3: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 2: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 1: match &= ~(*s1++ ^ *s2++); [[fallthrough]];
        case 0: default: ;
    }
    return is_all_zeros(match) ? match: popcnt512(match);
}
#endif
#endif

} // namespace detail

#ifdef _OPENMP
#define OMP_STATIC_SCHED32 _Pragma("omp parallel for schedule(static, 32)")
#else
#define OMP_STATIC_SCHED32
#endif

std::string path2cmd(const std::string &path) {
    if(std::equal(path.rbegin(), path.rbegin() + 3, "zg.")) return "gzip -dc "s + path;
    if(std::equal(path.rbegin(), path.rbegin() + 3, "zx.")) return "xz -dc "s + path;
    if(std::equal(path.rbegin(), path.rbegin() + 4, "2zb.")) return "bzip2 -dc "s + path;
    if(std::equal(path.rbegin(), path.rbegin() + 4, "tsz.")) return "zstd -dc "s + path;
    return ""s;
}

struct CompressedRet: public std::tuple<void *, long double, long double> {
    using super = std::tuple<void *, long double, long double>;
    std::unique_ptr<uint8_t[]> up;
    size_t nbytes = 0;
    bool ismapped = 0;
    CompressedRet &a(long double v) {
        std::get<1>(*this) = v;
        return *this;
    }
    CompressedRet &b(long double v) {
        std::get<2>(*this) = v;
        return *this;
    }
    long double a() const {return std::get<1>(*this);}
    long double b() const {return std::get<2>(*this);}
    CompressedRet(): super{nullptr, 0.L, 0.L} {
    }
    CompressedRet(CompressedRet &&o): super(static_cast<super>(o)), up(std::move(up)), nbytes(o.nbytes), ismapped(o.ismapped) {
        o.nbytes = 0; o.ismapped = 0; std::get<1>(o) = 0.; std::get<2>(o) = 0.;
        std::get<0>(o) = 0;
    }
    CompressedRet(const CompressedRet &o) = delete;
};

INLINE void *ptr_roundup(void *ptr) {
    uint8_t *tv = static_cast<uint8_t *>(ptr);
    if(auto rem = (uint64_t)tv & 63; rem) tv += (64 - rem);
    ptr = static_cast<void *>(tv);
    return ptr;
}

void make_compressed(CompressedRet &ret, int truncation_method, double fd, const mm::vector<RegT> &sigs, const mm::vector<uint64_t> &kmers, bool is_edit_distance, long double a=-1., long double b=-1., const bool sketch_compressed_set=0) {
    if(fd >= sizeof(RegT)) return;
    ret.nbytes = fd * sigs.size();
    if(fd == 0. && std::fmod(fd * sigs.size(), 1.)) THROW_EXCEPTION(std::runtime_error("Can't do nibble registers without an even number of signatures"));
    auto &compressed_reps = std::get<0>(ret);
    if(sketch_compressed_set) {
        compressed_reps = static_cast<void *>(const_cast<RegT *>(sigs.data()));
        ret.a(a);
        ret.b(b);
        assert(std::min(a, b) > 0.L);
        return;
    }
    ret.up.reset(new uint8_t[ret.nbytes + 63]);
    compressed_reps = ptr_roundup(static_cast<void *>(ret.up.get()));
    const size_t nsigs = sigs.size();
    assert(fd != 0.5 || sigs.size() % 2 == 0);
    if(is_edit_distance) {
        const uint64_t *sptr = (const uint64_t *)sigs.data();
        if(fd == 0.5) {
            uint8_t *ptr = (uint8_t *)compressed_reps;
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs / 2; ++i) {
                ptr[i] = (sptr[2 * i] & 0xfu) | ((sptr[2 * i + 1] & 0xfu) << 4);
            }
        } else {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs; ++i) {
                const auto sig = sptr[i];
                if(fd == 8) ((uint64_t *)compressed_reps)[i] = sig;
                else if(fd == 4) ((uint32_t *)compressed_reps)[i] = sig;
                else if(fd == 2) ((uint16_t *)compressed_reps)[i] = sig;
                else ((uint8_t *)compressed_reps)[i] = sig;
            }
        }
    } else if(truncation_method <= 0) {
        long double logbinv;
        const long double q = fd == 1. ? 254.3: fd == 2. ? 65534: fd == 4. ? 4294967294: fd == 8. ? 18446744073709551615. : fd == 0.5 ? 15.4: -1.;
        if(a <= 0. or b <= 0.) {
            RegT minreg = std::numeric_limits<RegT>::max(), maxreg = -std::numeric_limits<RegT>::max();
            OMP_ONLY(_Pragma("omp parallel for simd reduction(min:minreg) reduction(max:maxreg)"))
            for(size_t i = 0; i < nsigs; ++i) {
                const auto v = sigs[i];
                if(v <= 0 || v == std::numeric_limits<RegT>::max()) continue;
                minreg = std::min(minreg, v);
                maxreg = std::max(maxreg, v);
            }
            assert(q > 0.);
            const auto tailored = sketch::CSetSketch<RegT>::optimal_parameters(minreg, maxreg, q);
            b = tailored.first;
            a = tailored.second;
            std::fprintf(stderr, "Truncated via setsketch, a = %0.20Lg and b = %0.24Lg from min, max regs %Lg, %Lg\n", a, b, static_cast<long double>(minreg), static_cast<long double>(maxreg));
        }
        if(a == 0. || std::isinf(b)) {
            std::fprintf(stderr, "Note: setsketch compression yielded infinite value; falling back to b-bit compression\n");
            return make_compressed(ret, 1, fd, sigs, kmers, is_edit_distance);
        }
        ret.a(a); ret.b(b);
        logbinv = 1.L / std::log1p(b - 1.L);
        if(fd == 0.5) {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs / 2; ++i) {
                const auto lowerv = std::max((1.L - std::log(sigs[2 * i] / a) * logbinv), 0.L);
                const auto higherv = std::max((1.L - std::log(sigs[2 * i + 1] / a) * logbinv), 0.L);
                const uint8_t lower_half = std::max(0, std::min(int(q) + 1, static_cast<int>(lowerv)));
                const uint8_t upper_half = std::max(0, std::min(int(q) + 1, static_cast<int>(higherv))) << 4;
                ((uint8_t *)compressed_reps)[i] = lower_half | upper_half;
            }
        } else {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs; ++i) {
                const long double sub = 1.L - std::log(static_cast<long double>(sigs[i]) / a) * logbinv;
                if(fd == 8) {
                    ((uint64_t *)compressed_reps)[i] = std::min(uint64_t(q + 1), uint64_t(sub));
                } else {
                    const int64_t isub = std::max(int64_t(0), std::min(int64_t(q + 1), static_cast<int64_t>(sub)));
                    if(fd == 4)      ((uint32_t *)compressed_reps)[i] = isub;
                    else if(fd == 2) ((uint16_t *)compressed_reps)[i] = isub;
                    else             ((uint8_t *)compressed_reps)[i] = isub;
                }
            }
        }
    } else {
        //bbit:
        auto getsig = [kne=!kmers.empty(),&sigs,&kmers](auto x) -> uint64_t {
            return kne ? sketch::hash::WangHash::hash(kmers[x]): reg2sig(sigs[x]);
        };
        if(fd == 0.5) {
            uint8_t *ptr = (uint8_t *)compressed_reps;
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs / 2; ++i) {
                const uint64_t sig1 = getsig(2 * i), sig2 = getsig(2 * i + 1);
                ptr[i] = (sig1 & 0xfu) | ((sig2 & 0xfu) << 4);
            }
        } else {
            static constexpr int shift [] {0, 58, 48, 0, 32, 0, 0, 0, 0, 0};
            static_assert(shift[1] == 58, "Shift 1 must be 58");
            static_assert(shift[2] == 48, "Shift 2 must be 48");
            static_assert(shift[4] == 32, "Shift 4 must be 32");
            static_assert(shift[8] == 0, "Shift 8 must be 0");
            const int myshift = shift[int(fd)];
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs; ++i) {
                const uint64_t sig = getsig(i) >> myshift;
                if(fd == 8) ((uint64_t *)compressed_reps)[i] = sig;
                else if(fd == 4) ((uint32_t *)compressed_reps)[i] = sig;
                else if(fd == 2) ((uint16_t *)compressed_reps)[i] = sig;
                else ((uint8_t *)compressed_reps)[i] = sig;
            }
        }
    }
}
static inline long double g_b(long double b, long double arg) {
    return (1.L - std::pow(b, -arg)) / (1.L - 1.L / b);
}

LSHDistType compare(const Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j) {
    long double ret = std::numeric_limits<LSHDistType>::max();
    const long double lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
    const long double invdenom = 1.L / opts.sketchsize_;
    auto sim2dist = [poisson_mult=-1. / std::max(1, opts.k_)](auto x) -> double {if(x) return std::log(2. * x / (1. + x)) * poisson_mult; return std::numeric_limits<double>::infinity();};
    if(opts.compressed_ptr_) {
        const bool bbit_c = opts.truncation_method_ > 0;
        std::pair<uint64_t, uint64_t> res{0, 0};
        if(bbit_c) {
            auto &equal_regs = std::get<0>(res);
            switch(int(2. * opts.fd_level_)) {

#define CASEPOW2 \
                CASE_ENTRY(16, uint64_t)\
                CASE_ENTRY(8, uint32_t)\
                CASE_ENTRY(4, uint16_t)\
                CASE_ENTRY(2, uint8_t)
#define CASE_ENTRY(v, TYPE)\
case v: {TYPE *ptr = static_cast<TYPE *>(opts.compressed_ptr_); equal_regs = sketch::eq::count_eq(ptr + i * opts.sketchsize_, ptr + j * opts.sketchsize_, opts.sketchsize_);} break;
                CASEPOW2
#undef CASE_ENTRY
                case 1: {
                    uint8_t *ptr = static_cast<uint8_t *>(opts.compressed_ptr_);
                    equal_regs = sketch::eq::count_eq_nibbles(ptr + i * opts.sketchsize_ / 2, ptr + j * opts.sketchsize_ / 2, opts.sketchsize_);
                    break;
                }
                default: __builtin_unreachable();
            }
            //std::fprintf(stderr, "%zu equal registers out of %zu\n", size_t(equal_regs), opts.sketchsize_);
        } else {
            switch(int(2. * opts.fd_level_)) {
#define CASE_ENTRY(v, TYPE)\
case v: {\
    TYPE *ptr = static_cast<TYPE *>(opts.compressed_ptr_);\
    res = sketch::eq::count_gtlt(ptr + i * opts.sketchsize_, ptr + j * opts.sketchsize_, opts.sketchsize_);\
    } break;
                CASEPOW2
#undef CASE_ENTRY
#undef CASEPOW2
                case 1: {
                    uint8_t *ptr = static_cast<uint8_t *>(opts.compressed_ptr_);
                    res = sketch::eq::count_gtlt_nibbles(ptr + i * opts.sketchsize_ / 2, ptr + j * opts.sketchsize_ / 2, opts.sketchsize_);
                    break;
                }
                default: __builtin_unreachable();
            }
        }
        if(bbit_c) {
            // ret = ((num / denom) - (1. / 2^b)) / (1. - 1. / 2^b);
            // maps equality to 1 and down-estimates for account for collisions
            const long double b2pow = -std::ldexp(1.L, -static_cast<int>(opts.fd_level_ * 8.));
            ret = std::max(0.L, std::fma(res.first, invdenom, b2pow) / (1.L + b2pow));
            if(opts.measure_ == INTERSECTION)
                ret *= std::max((lhcard + rhcard) / (2.L - (1.L - ret)), 0.L);
            else if(opts.measure_ == CONTAINMENT)
                ret = std::max((lhcard + rhcard) / (2.L - (1.L - ret)), 0.L) * ret / lhcard;
            else if(opts.measure_ == POISSON_LLR)
                ret = sim2dist(ret);
            else if(opts.measure_ == SYMMETRIC_CONTAINMENT)
                ret = std::max((lhcard + rhcard) / (2.L - (1.L - ret)), 0.L) * ret / std::min(lhcard, rhcard);
        } else {
            long double alpha = res.first * invdenom;
            long double beta = res.second * invdenom;
            const long double b = opts.compressed_b_;
            static_assert(sizeof(b) == 16, "b must be a long double");
            long double mu;
            if(opts.fd_level_ < sizeof(RegT)) {
                alpha = g_b(b, alpha);
                beta = g_b(b, beta);
            }
            //std::fprintf(stderr, "alpha gt: %zu. beta gt: %zu, Alpha %Lg, Beta %Lg\n", res.first, res.second, alpha, beta);
            VERBOSE_ONLY(std::fprintf(stderr, "Alpha: %Lg. Beta: %Lg\n", alpha, beta);)
            if(alpha + beta >= 1.) {
                mu = lhcard + rhcard;
            } else {
                mu = std::max((lhcard + rhcard) / (2.L - alpha - beta), 0.L);
            }
            auto triple = std::make_tuple(alpha, beta, mu);
            //std::fprintf(stderr, "Triple: %Lg apha, %Lg beta, %Lg mu\n", std::get<0>(triple), std::get<1>(triple), std::get<2>(triple));
            ret = std::max(1.L - (std::get<0>(triple) + std::get<1>(triple)), 0.L);
            if(opts.measure_ == INTERSECTION)
                ret *= mu;
            else if(opts.measure_ == CONTAINMENT)
                ret = (ret * mu) / lhcard;
            else if(opts.measure_ == SYMMETRIC_CONTAINMENT)
                ret = (ret * mu) / std::min(lhcard, rhcard);
            else if(opts.measure_ == POISSON_LLR)
                ret = sim2dist(ret);
        }
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE && (opts.exact_kmer_dist_ || opts.measure_ == M_EDIT_DISTANCE)) {
        //std::fprintf(stderr, "Return exact distance for edit distance between sequences\n");
        assert(result.sequences_.size() > std::max(i, j) || !std::fprintf(stderr, "Expected sequences to be non-null for exact edit distance calculation (%zu vs %zu/%zu)\n", result.sequences_.size(), i, j));
        return levenshteinSSE::levenshtein(result.sequences_[i], result.sequences_[j]);
    } else if(opts.kmer_result_ <= FULL_SETSKETCH) {
        const RegT *lhsrc = &result.signatures_[opts.sketchsize_ * i], *rhsrc = &result.signatures_[opts.sketchsize_ * j];
        if(opts.sspace_ == SPACE_SET && opts.truncation_method_ <= 0) {
            auto gtlt = sketch::eq::count_gtlt(lhsrc, rhsrc, opts.sketchsize_);
            long double alpha, beta, eq, lhcard, ucard, rhcard;
            alpha = gtlt.first * invdenom;
            beta = gtlt.second * invdenom;
            lhcard = result.cardinalities_[i], rhcard = result.cardinalities_[j];
            eq = (1. - alpha - beta);
            if(eq <= 0.)
                return opts.measure_ != POISSON_LLR ? 0.: std::numeric_limits<double>::max();
            ucard = std::max(lhcard + rhcard / (2.L - alpha - beta), 0.L);
            LSHDistType isz = ucard * eq, sim = eq;
            ret = opts.measure_ == SIMILARITY ? sim
                : opts.measure_ == INTERSECTION ? isz
                : opts.measure_ == SYMMETRIC_CONTAINMENT ? isz / (std::min(lhcard, rhcard))
                : opts.measure_ == POISSON_LLR ? sim2dist(sim): LSHDistType(-1);
            assert(ret >= 0. || !std::fprintf(stderr, "measure: %s. sim: %g. isz: %g\n", to_string(opts.measure_).data(), sim, isz));
        } else {
            const RegT *sptr = result.signatures_.data();
            if constexpr(sizeof(RegT) == 8) {
                // If RegT are the same size as k-mers, compare the k-mers themselves
                // instead of the doubles
                // Since we're only comparing for equality, this can only improve accuracy
                if(result.kmers_.size() == result.signatures_.size() && !opts.use128()) {
                    DBG_ONLY(std::fprintf(stderr, "Comparing k-mers sampled rather than the items themselves. This should be more specific, since there is 0 chance of collisions.\n");)
                    sptr = reinterpret_cast<const RegT *>(result.kmers_.data());
                }
            }
            const auto neq = sketch::eq::count_eq(&sptr[opts.sketchsize_ * i], &sptr[opts.sketchsize_ * j], opts.sketchsize_);
            ret = invdenom * neq;
            if(opts.measure_ == INTERSECTION) {
                ret *= std::max((lhcard + rhcard) / (1.L + ret), 0.L);
            } else if(opts.measure_ == SYMMETRIC_CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1.L + ret), 0.L) / std::min(lhcard, rhcard);
            else if(opts.measure_ == CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1.L + ret), 0.L) / lhcard;
            else if(opts.measure_ == POISSON_LLR) ret = sim2dist(ret);
        }
    } else {
#define CORRECT_RES(res, measure, lhc, rhc)\
            if(measure == SYMMETRIC_CONTAINMENT) \
                res = res / std::min(lhc, rhc);\
            else if(measure == POISSON_LLR || measure == SIMILARITY)\
                res = res / (lhc + rhc - res);\
            else if(measure == CONTAINMENT) res /= lhc;\
            ret = res;
        const std::string &lpath = result.destination_files_[i], &rpath = result.destination_files_[j];
        if(lpath.empty() || rpath.empty()) THROW_EXCEPTION(std::runtime_error("Destination files for k-mers empty -- cannot load from disk"));
        std::FILE *lhk = 0, *rhk = 0, *lhn = 0, *rhn = 0;
        std::string lcmd = path2cmd(lpath), rcmd = path2cmd(rpath), lkcmd, rkcmd;
        lhk = lcmd.empty() ? bfopen(lpath.data(), "rb"): ::popen(lcmd.data(), "r");
        rhk = rcmd.empty() ? bfopen(rpath.data(), "rb"): ::popen(rcmd.data(), "r");
        if(lhk == 0) THROW_EXCEPTION(std::runtime_error("Failed to "s + (lcmd.empty() ? " read from "s + lpath: "Run cmd '"s + lcmd + "'")));
        if(rhk == 0) THROW_EXCEPTION(std::runtime_error("Failed to "s + (rcmd.empty() ? " read from "s + rpath: "Run cmd '"s + rcmd + "'")));
        if(result.kmercountfiles_.size()) {
            lkcmd = path2cmd(result.kmercountfiles_[i]);
            rkcmd = path2cmd(result.kmercountfiles_[j]);
            lhn = lkcmd.empty() ? bfopen(result.kmercountfiles_[i].data(), "rb"): ::popen(lkcmd.data(), "r");
            if(rkcmd.empty()) {
                rhn = bfopen(result.kmercountfiles_[j].data(), "rb");
            } else {
                rhn = ::popen(rkcmd.data(), "r");
            }
            if(lhn == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to read from "s + result.kmercountfiles_[i]));
            if(rhn == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to read from "s + result.kmercountfiles_[j]));
        }
        if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
            if(opts.exact_kmer_dist_) {
                auto [edit_dist, max_edit_dist] = mmer_edit_distance(lhk, rhk, opts.use128());
                ret = opts.measure_ == M_EDIT_DISTANCE ? edit_dist: max_edit_dist - edit_dist;
            } else {
                ret = hamming_compare_f64(lhk, rhk);
            }
        } else {
            double lhc, rhc;
            std::fread(&lhc, sizeof(lhc), 1, lhk);
            std::fread(&rhc, sizeof(rhc), 1, rhk);
            const auto [isz_size, union_size] = weighted_compare(lhk, rhk, lhn, rhn, lhc, rhc, opts.use128());
            double res = isz_size;
            CORRECT_RES(res, opts.measure_, lhc, rhc)
        }
        if(lcmd.empty()) std::fclose(lhk); else ::pclose(lhk);
        if(rcmd.empty()) std::fclose(rhk); else ::pclose(rhk);
        if(lhn) {
            if(lkcmd.empty()) std::fclose(lhn); else ::pclose(lhn);
        }
        if(rhn) {
            if(rkcmd.empty()) std::fclose(rhn); else ::pclose(rhn);
        }
#undef CORRECT_RES
        // Compare exact representations, not compressed shrunk
    }
    if(std::isnan(ret) || std::isinf(ret)) ret = std::numeric_limits<decltype(ret)>::max();
    return ret;
}

template<typename MHT>
inline size_t densify(MHT *minhashes, uint64_t *const kmers, const size_t sketchsize, const schism::Schismatic<uint64_t> &div, const MHT empty=MHT(0))
{
    const long long unsigned int ne = simdcount(minhashes, sketchsize, empty);
    if(ne == sketchsize  || ne == 0) return ne;
    for(size_t i = 0; i < sketchsize; ++i) {
        if(minhashes[i] != empty) continue;
        uint64_t j = i;
        uint64_t rng = i;
        while(minhashes[j] == empty) {
            static constexpr uint32_t PRIMEMOD = 4294967291u;
            auto a = (wy::wyhash64_stateless(&rng) % PRIMEMOD), b = (wy::wyhash64_stateless(&rng) % PRIMEMOD), c = (wy::wyhash64_stateless(&rng) % PRIMEMOD);
            j = div.mod((((a * b) % PRIMEMOD + c) % PRIMEMOD));
        }
        minhashes[i] = minhashes[j];
        if(kmers) kmers[i] = kmers[j];
    }
    assert(std::find(minhashes, minhashes + sketchsize, empty) == minhashes + sketchsize);
    return ne;
}

void cmp_core(const Dashing2DistOptions &opts, SketchingResult &result) {
    if(opts.sketch_compressed() && opts.truncate_mode() != 0) {
        THROW_EXCEPTION(std::invalid_argument("Can't use truncated setsketch generation with bbit signatures. Omit --bbit-sigs or --setsketch-ab"));
    }
    if(opts.sketch_compressed() && opts.save_kmers()) {
        THROW_EXCEPTION(std::invalid_argument("Can't use truncated setsketch generation --save-kmers. Omit --save-kmers or --setsketch-ab"));
    }
    // We handle some details before dispatching the final comparison code
    // First, we compute cardinalities for sets/multisets
    // and then we densify one-permutation minhashing
    if(opts.kmer_result_ >= FULL_MMER_SET && result.cardinalities_.size() && result.cardinalities_.front() < 0.) {
        const size_t n = result.cardinalities_.size();
        if(opts.parse_by_seq_)
            THROW_EXCEPTION(std::runtime_error("Not yet supported"));
        OMP_PFOR
        for(size_t i = 0; i < n; ++i) {
            int ft;
            std::FILE *ifp = nullptr;
            std::string fn = opts.kmer_result_ == FULL_MMER_SET ? result.kmerfiles_.at(i): result.destination_files_.at(i);
            std::string cmd;
            if(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                if(!check_compressed(fn, ft)) throw std::runtime_error(std::string("Missing kmerfile or destination file: ") + fn);
                if(endswith(fn, ".xz")) ft = 2;
                else if(endswith(fn, ".gz")) ft = 1;
                else ft = 0;
                cmd = std::string(ft == 0 ? "": ft == 1 ? "gzip -dc ": ft == 2 ? "xz -dc " : "unknowncommand");
                if(cmd.empty()) {
                    struct stat st;
                    if(::stat(fn.data(), &st)) {
                        perror((std::string("Failed to stat ") + fn).data());
                        std::exit(1);
                    }
                    result.cardinalities_[i] = st.st_size / sizeof(uint64_t);
                    assert(st.st_size % sizeof(uint64_t) == 0);
                } else {
                    cmd = cmd + " " + fn;
                    ifp = ::popen(cmd.data(), "r");
                    if(!ifp)
                        THROW_EXCEPTION(std::runtime_error("Failed to read from '"s + cmd + "' for file " + fn));
                    size_t c = 0;
                    for(uint64_t x;std::fread(&x, sizeof(x), 1, ifp) == 1u;++c)
                    result.cardinalities_[i] = c;
                }
            } else if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                if(!check_compressed(result.kmercountfiles_[i], ft)) throw std::runtime_error("Missing kmercountfile");
                if(endswith(result.kmercountfiles_[i], ".xz")) ft = 2;
                else if(endswith(result.kmercountfiles_[i], ".gz")) ft = 1;
                else ft = 0;
                cmd = ft == 0 ? ""s: ft == 1 ? "gzip -dc "s: ft == 2 ? "xz -dc "s : "unknowncommand"s;
                if(cmd == "unknowncommand") THROW_EXCEPTION(std::runtime_error("Failure"));
                if(cmd.empty()) {
                    ifp = bfopen(fn.data(), "rb");
                } else {
                    cmd += result.kmercountfiles_[i];
                    ifp = ::popen(cmd.data(), "r");
                }
                if(0 == ifp)
                    THROW_EXCEPTION(std::runtime_error(cmd.empty() ? "Parsing "s + fn : ": Command '"s + cmd + "' failed."));
                double x, c, s;
                for(x = c = s = 0.;std::fread(&x, sizeof(x), 1, ifp) == 1u;sketch::kahan::update(s, c, x));
                result.cardinalities_[i] = s;
            }
            if(ifp) {
                if(cmd.empty()) std::fclose(ifp);
                else ::pclose(ifp);
            }
        }
    }
    if(opts.kmer_result_ == ONE_PERM) {
        const schism::Schismatic<uint64_t> sd(opts.sketchsize_);
        const size_t n = result.cardinalities_.size();
        uint64_t *const kp = result.kmers_.size() ? result.kmers_.data(): (uint64_t *)nullptr;
        size_t totaldens = 0;
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 32) reduction(+:totaldens)
#endif
        for(size_t i = 0; i < n;++i) {
            totaldens += densify(&result.signatures_[opts.sketchsize_ * i],
                                 kp ? &kp[opts.sketchsize_ * i]: kp,
                                 opts.sketchsize_, sd);
        }
        if(totaldens > 0) std::fprintf(stderr, "Densified a total of %zu/%zu entries\n", totaldens, opts.sketchsize_ * n);
    }
    //std::fprintf(stderr, "Handled generating needed cardinalities\n");
    // Calculate cardinalities if they haven't been
    VERBOSE_ONLY(
    std::fprintf(stderr, "Beginning cmp_core with options: \n");
        if(opts.sspace_ == SPACE_SET) {
            std::fprintf(stderr, "Comparing sets\n");
        } else if(opts.sspace_ == SPACE_MULTISET) {
            std::fprintf(stderr, "Comparing multisets\n");
        } else if(opts.sspace_ == SPACE_PSET) {
            std::fprintf(stderr, "Comparing discrete probability distributions\n");
        } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            std::fprintf(stderr, "Comparing items in edit distance space\n");
        }
        std::fprintf(stderr, "Result type: %s\n", to_string(opts.kmer_result_).data());
    )
    if(opts.kmer_result_ <= FULL_MMER_SET && opts.fd_level_ < sizeof(RegT)) {
        if(result.signatures_.empty()) THROW_EXCEPTION(std::runtime_error("Empty signatures; trying to compress registers but don't have any"));
    }
    CompressedRet cret;
    make_compressed(cret, opts.truncation_method_, opts.fd_level_, result.signatures_, result.kmers_, opts.sspace_ == SPACE_EDIT_DISTANCE, opts.compressed_a_, opts.compressed_b_, opts.sketch_compressed_set);
    std::tie(opts.compressed_ptr_, opts.compressed_a_, opts.compressed_b_) = cret;
    if(opts.output_kind_ <= ASYMMETRIC_ALL_PAIRS || opts.output_kind_ == PANEL) {
        emit_rectangular(opts, result);
        return;
    }
    // This is LSH-index assisted KNN graphs +
    // thresholded nn graphs

    // Step 1: Build LSH Index
    std::vector<uint64_t> nperhashes;
    while(nperhashes.size() < opts.nLSH) {
        nperhashes.emplace_back(nperhashes.size() < 3 ? (1ull << nperhashes.size()): static_cast<unsigned long long>(nperhashes.size() * 2));
    }
    std::vector<uint64_t> nperrows(nperhashes.size());
    for(size_t i = 0; i < nperhashes.size(); ++i) {
        const auto nh = nperhashes[i];
        auto &np = nperrows[i];
        if(nh < 2) {
            np = opts.sketchsize_ / nh;
        } else if(nh <= 4) {
            np = opts.sketchsize_ * 8 / nh;
        } else if(nh <= 6) {
            np = opts.sketchsize_  * 16 / nh;
        } else np = opts.sketchsize_ * 32 / nh;
    }
    using SSI = SetSketchIndex<LSHIDType, LSHIDType>;
    SSI idx(opts.kmer_result_ < FULL_MMER_SET ? SSI(opts.sketchsize_, nperhashes, nperrows): SSI());


    // Step 2: Build nearest-neighbor candidate table
    if(opts.output_kind_ == KNN_GRAPH || opts.output_kind_ == NN_GRAPH_THRESHOLD) {
        const bool exact_knn = std::getenv("EXACT_KNN");
        std::vector<pqueue> neighbor_lists = exact_knn ? build_exact_graph(idx, opts, result): build_index(idx, opts, result);
        if(!exact_knn) {
            refine_results(neighbor_lists, opts, result);
        } else if(!distance(opts.measure_)) {
            OMP_PFOR
            for(size_t i = 0; i < neighbor_lists.size(); ++i) {
                auto &nl = neighbor_lists[i];
                auto beg = nl.begin(), e = nl.end();
                std::transform(beg, e, beg, [&](PairT x) {return PairT{-x.first, x.second};});
            }
        }
        emit_neighbors(neighbor_lists, opts, result);
    } else if(opts.output_kind_ == DEDUP) {
        // The ID corresponds to the representative of a cluster;
        // Constituents is a vector of IDs per cluster;
        auto [ids, constituents] = dedup_core(idx, opts, result);
        dedup_emit(ids, constituents, opts, result);
    }
}

} // namespace dashing2
