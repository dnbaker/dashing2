#include "cmp_main.h"
#include "sketch/count_eq.h"
//#include "sketch/ssi.h"
#include "index_build.h"
#include "refine.h"


namespace dashing2 {
//using sketch::lsh::SetSketchIndex;
static INLINE uint64_t reg2sig(RegT x) {
    uint64_t seed = 0;
    CONST_IF(sizeof(RegT) <= 8) {
        std::memcpy(&seed, &x, sizeof(seed));
        return wy::wyhash64_stateless(&seed);
    } else {
        std::memcpy(&seed, &x, sizeof(seed));
        uint64_t nextseed = wy::wyhash64_stateless(&seed);
        nextseed ^= ((uint64_t *)&x)[1];
        return wy::wyhash64_stateless(&nextseed);
    }
}

std::tuple<void *, double, double>
make_compressed(int truncation_method, double fd, const std::vector<RegT> &sigs) {
    std::tuple<void *, double, double> ret = {nullptr, 0., 0.};
    auto &compressed_reps = std::get<0>(ret);
    if(fd >= sizeof(RegT)) return ret;
    size_t mem = fd * sigs.size();
    if(fd == 0. && std::fmod(fd * sigs.size(), 1.)) throw std::runtime_error("Can't do nibble registers without an even number of signatures");
    if(posix_memalign((void **)&compressed_reps, 64, mem)) throw std::bad_alloc();
    const size_t nsigs = sigs.size();
    if(fd == 0.5) {
        assert(sigs.size() % 2 == 0);
    }
    if(truncation_method > 0) {
        std::fprintf(stderr, "Performing %d-bit compression\n", int(fd * 8.));
        if(fd == 8) {
            std::transform(sigs.begin(), sigs.end(), (uint64_t *)compressed_reps, reg2sig);
        } else if(fd == 4) {
            std::transform(sigs.begin(), sigs.end(), (uint32_t *)compressed_reps, reg2sig);
        } else if(fd == 2) {
            std::transform(sigs.begin(), sigs.end(), (uint16_t *)compressed_reps, reg2sig);
        } else if(fd == 1) {
            std::transform(sigs.begin(), sigs.end(), (uint8_t *)compressed_reps, reg2sig);
        } else if(fd == 0.5) {
            uint8_t *ptr = (uint8_t *)compressed_reps;
            for(size_t i = 0; i < nsigs / 2; ++i) {
                auto sig1 = sigs[2 * i], sig2 = sigs[2 * i + 1];
                ptr[i] = (reg2sig(sig1) & 0xfu) | ((reg2sig(sig2) & 0xfu) << 4);
                //std::fprintf(stderr, "Double reg = %u\n", ptr[i], (reg2sig(sig1) & 0xfu), (reg2sig(sig2) & 0xfu) << 4);
            }
        } else __builtin_unreachable();
    } else {
        std::fprintf(stderr, "Performing logarithmic %d-bit compression\n", int(fd * 8.));
        RegT minreg = std::numeric_limits<RegT>::max(), maxreg = std::numeric_limits<RegT>::min();
        #pragma omp simd
        for(size_t i = 0; i < nsigs; ++i) {
            minreg = std::min(minreg, sigs[i]);
            maxreg = std::max(maxreg, sigs[i]);
        }
        //std::fprintf(stderr, "fd: %g\n", fd);
        double q = fd == 1. ? 254.3: fd == 2. ? 65534.3: fd == 4. ? 4294967294.3: fd == 8. ? 18446744073709551615. : fd == 0.5 ? 15.4: -1.;
        if(q < 0) throw 2;
        auto [b, a] = sketch::CSetSketch<RegT>::optimal_parameters(minreg, maxreg, q);
        std::fprintf(stderr, "Truncated via setsketch, a = %Lg and b = %Lg\n", a, b);
        std::get<1>(ret) = a;
        std::get<2>(ret) = b;
        const RegT logbinv = 1. / std::log1p(b - 1.);
        if(fd == 4) {
            std::transform(sigs.begin(), sigs.end(), (uint32_t *)compressed_reps,
                           [q,a,logbinv](RegT x) {
                 return std::max(uint32_t(0), std::min(uint32_t(q) + 1, static_cast<uint32_t>((1. - std::log(x / a) * logbinv))));
            });
        } else if(fd == 8) {
            std::transform(sigs.begin(), sigs.end(), (uint64_t *)compressed_reps,
                           [q,a,logbinv](RegT x)
            {
                 return std::max(int64_t(0), std::min(int64_t(q) + 1, static_cast<int64_t>((1. - std::log(x / a) * logbinv))));
            });
        } else if(fd == 2) {
            std::transform(sigs.begin(), sigs.end(), (uint16_t *)compressed_reps,
                           [q,a,logbinv](RegT x) {
                 return std::max(int64_t(0), std::min(int64_t(q) + 1, static_cast<int64_t>((1. - std::log(x / a) * logbinv))));
            });
        } else if(fd == 1) {
            std::transform(sigs.begin(), sigs.end(), (uint8_t *)compressed_reps,
                           [q,a,logbinv](RegT x) {
                 return std::max(int(0), std::min(int(q) + 1, static_cast<int>((1. - std::log(x / a) * logbinv))));
            });
        } else if(fd == 0.5) {
            uint8_t *ptr = (uint8_t *)compressed_reps;
            for(size_t i = 0; i < nsigs / 2; ++i) {
                uint8_t lower_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i] / a) * logbinv))));
                uint8_t upper_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i + 1] / a) * logbinv)))) << 4;
                ptr[i] = lower_half | upper_half;
            }
        }
    }
    return ret;
}
LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j) {
    LSHDistType ret = std::numeric_limits<LSHDistType>::max();
    const LSHDistType b2pow = -std::ldexp(1., -static_cast<int>(opts.fd_level_ * 8.));
    const LSHDistType ib2pow = 1. / (1. + b2pow);
    const LSHDistType invdenom = 1. / opts.sketchsize_;
    const LSHDistType lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
    const LSHDistType poisson_mult = 1. / std::max(1, opts.k_);
    if(opts.compressed_ptr_) {
        const bool bbit_c = opts.truncation_method_ > 0;
        std::pair<uint64_t, uint64_t> res{0, 0};
        if(bbit_c) {
            auto &equal_regs = std::get<0>(res);
            switch(int(2. * opts.fd_level_)) {

#define CASE_ENTRY(v, TYPE)\
case v: {TYPE *ptr = static_cast<TYPE *>(opts.compressed_ptr_); equal_regs = sketch::eq::count_eq(ptr + i * opts.sketchsize_, ptr + j * opts.sketchsize_, opts.sketchsize_);} break;
                CASE_ENTRY(16, uint64_t)
                CASE_ENTRY(8, uint32_t)
                CASE_ENTRY(4, uint16_t)
                CASE_ENTRY(2, uint8_t)
#undef CASE_ENTRY
                case 1: {
                    uint8_t *ptr = static_cast<uint8_t *>(opts.compressed_ptr_);
                    equal_regs = sketch::eq::count_eq_nibbles(ptr + i * opts.sketchsize_ / 2, ptr + j * opts.sketchsize_ / 2, opts.sketchsize_);
                    break;
                }
                default: __builtin_unreachable();
            }
        } else {
            switch(int(2. * opts.fd_level_)) {
#define CASE_ENTRY(v, TYPE)\
case v: {TYPE *ptr = static_cast<TYPE *>(opts.compressed_ptr_); res = sketch::eq::count_gtlt(ptr + i * opts.sketchsize_, ptr + j * opts.sketchsize_, opts.sketchsize_);} break;
                CASE_ENTRY(16, uint64_t)
                CASE_ENTRY(8, uint32_t)
                CASE_ENTRY(4, uint16_t)
                CASE_ENTRY(2, uint8_t)
#undef CASE_ENTRY
                case 1: {
                    uint8_t *ptr = static_cast<uint8_t *>(opts.compressed_ptr_);
                    res = sketch::eq::count_gtlt_nibbles(ptr + i * opts.sketchsize_ / 2, ptr + j * opts.sketchsize_ / 2, opts.sketchsize_);
                    std::fprintf(stderr, "gt/lt/eq: %zu/%zu/%zu\n", size_t(res.first), size_t(res.second), size_t(opts.sketchsize_ - res.first - res.second));
                    break;
                }
                default: __builtin_unreachable();
            }
        }
        if(bbit_c) {
            // ret = ((num / denom) - (1. / 2^b)) / (1. - 1. / 2^b);
            // maps equality to 1 and down-estimates for account for collisions
            ret = std::max(0., std::fma(res.first, invdenom, b2pow) * ib2pow);
            if(opts.measure_ == INTERSECTION)
                ret *= std::max((lhcard + rhcard) / (2. - (1. - ret)), 0.);
            else if(opts.measure_ == CONTAINMENT)
                ret = std::max((lhcard + rhcard) / (2. - (1. - ret)), 0.) * ret / lhcard;
            else if(opts.measure_ == POISSON_LLR)
                ret = -std::log(ret) * poisson_mult;
            else if(opts.measure_ == SYMMETRIC_CONTAINMENT)
                ret = std::max((lhcard + rhcard) / (2. - (1. - ret)), 0.) * ret / std::min(lhcard, rhcard);
        } else {
            auto alpha = res.first * invdenom;
            auto beta = res.second * invdenom;
            const LSHDistType b = opts.compressed_b_;
            LSHDistType mu;
            if(opts.fd_level_ < sizeof(RegT)) {
                alpha = sketch::setsketch::g_b(b, alpha);
                beta = sketch::setsketch::g_b(b, beta);
            }
            if(alpha + beta >= 1.) {
                mu = lhcard + rhcard;
            } else {
                mu = std::max((lhcard + rhcard) / (2. - alpha - beta), 0.);
            }
            auto triple = std::make_tuple(alpha, beta, mu);
            ret = std::max(1. - (std::get<0>(triple) + std::get<1>(triple)), 0.);
            if(opts.measure_ == INTERSECTION)
                ret *= mu;
            else if(opts.measure_ == CONTAINMENT)
                ret = (ret * mu) / lhcard;
            else if(opts.measure_ == SYMMETRIC_CONTAINMENT)
                ret = (ret * mu) / std::min(lhcard, rhcard);
            else if(opts.measure_ == POISSON_LLR)
                ret = -std::log(ret) * poisson_mult;
        }
    } else if(opts.kmer_result_ <= FULL_SETSKETCH) {
        const RegT *lhsrc = &result.signatures_[opts.sketchsize_ * i], *rhsrc = &result.signatures_[opts.sketchsize_ * j];
        if(opts.sspace_ == SPACE_SET) {
            auto gtlt = sketch::eq::count_gtlt(lhsrc, rhsrc, opts.sketchsize_);
            LSHDistType alpha = gtlt.first * invdenom;
            LSHDistType beta = gtlt.second * invdenom;
            LSHDistType eq = (1. - alpha - beta);
            LSHDistType lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
            LSHDistType ucard = alpha + beta >= 1. ? lhcard + rhcard: std::max((lhcard + rhcard) / (2. - alpha - beta), 0.);
            //std::fprintf(stderr, "alpha %g, beta %g, eq = %g, lhcard %g, rhcard %g, ucard %g\n", alpha, beta, eq, lhcard, rhcard, ucard);
            //std::fprintf(stderr, "gtlt %zu/%zu\n", gtlt.first, gtlt.second);
            LSHDistType isz = ucard * eq;
            LSHDistType sim = isz / ucard;
            ret = opts.measure_ == SIMILARITY ? sim
                : opts.measure_ == INTERSECTION ? isz
                : opts.measure_ == SYMMETRIC_CONTAINMENT ? isz / (std::min(lhcard, rhcard))
                : opts.measure_ == POISSON_LLR ? -std::log(sim) * poisson_mult: -1.;
            assert(ret >= 0.);
        } else {
            ret = invdenom * sketch::eq::count_eq(&result.signatures_[opts.sketchsize_ * i], &result.signatures_[opts.sketchsize_ * j], opts.sketchsize_);
            if(opts.measure_ == INTERSECTION) ret *= std::max((lhcard + rhcard) / (1. + ret), 0.);
            else if(opts.measure_ == SYMMETRIC_CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1. + ret), 0.) / std::min(lhcard, rhcard);
            else if(opts.measure_ == CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1. + ret), 0.) / lhcard;
            else if(opts.measure_ == POISSON_LLR) ret = -std::log(ret) * poisson_mult;
        }
    } else {
        throw std::runtime_error("Not yet implemented: pairwise distances for full hash sets or count dictionaries; only sketching comparisons is allowed currently");
        // Compare exact representations, not compressed shrunk
    }
    if(std::isnan(ret) || std::isinf(ret)) ret = std::numeric_limits<double>::max();
    return ret;
}
void emit_all_pairs(Dashing2DistOptions &opts, const SketchingResult &result) {
    const size_t ns = result.names_.size();
#ifndef NDEBUG
    for(size_t i = 0; i < ns; ++i) {
        std::fprintf(stderr, "name %s for index %zu\n", result.names_[i].data(), i);
    }
#endif
    std::FILE *ofp = opts.outfile_path_.empty() ? stdout: std::fopen(opts.outfile_path_.data(), "wb");
    if(!ofp) throw std::runtime_error(std::string("Failed to open path at ") + opts.outfile_path_);
    std::setvbuf(ofp, nullptr, _IOFBF, 1<<17);
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::thread sub;
    for(size_t i = 0; i < ns; ++i) {
        // TODO: batch queries together for cache efficiency (distmat::parallel_fill for an example)
        size_t nelem = asym ? ns: ns - i - 1;
        std::unique_ptr<float[]> dat(new float[nelem]);
        const size_t diff = asym ? 0: i + 1;
        const size_t nwritten = asym ? ns: ns - i - 1;
        OMP_PFOR_DYN
        for(size_t start = asym ? 0: i + 1;start < ns; ++start) {
            auto v = compare(opts, result, i, start);
            //std::fprintf(stderr, "Calling %zu, %zu (i < ns = %zu)\n", i, start, ns);
            dat[start - diff] = v;
        }
        if(sub.joinable()) sub.join();
        // TODO: turn the unique pointer to an output queue
        sub = std::thread([dat=std::move(dat),nwritten,ofp,i]() {
            if(std::fwrite(dat.get(), sizeof(float), nwritten, ofp) != nwritten)
                throw std::runtime_error(std::string("Failed to write row ") + std::to_string(i) + " to disk");
        }); 
    }
    if(ofp != stdout) std::fclose(ofp);
    sub.join();
}
void cmp_core(Dashing2DistOptions &opts, const SketchingResult &result) {
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
    if(opts.kmer_result_ <= FULL_MMER_SET && opts.fd_level_ < sizeof(RegT)) {
        if(result.signatures_.empty()) throw std::runtime_error("Empty signatures; trying to compress registers but don't have any");
    }
    std::fprintf(stderr, "Result type: %s\n", to_string(opts.kmer_result_).data());
    std::tie(opts.compressed_ptr_, opts.compressed_a_, opts.compressed_b_) = make_compressed(opts.truncation_method_, opts.fd_level_, result.signatures_);
    if(opts.output_kind_ <= ASYMMETRIC_ALL_PAIRS) {
        emit_all_pairs(opts, result);
    } else {
        // Step 1: Build LSH Index
        std::vector<uint64_t> nperhashes{1, 2, 3, 4, 6, 8, 16};
        std::vector<uint64_t> nperrows(nperhashes.size());
        for(size_t i = 0; i < nperhashes.size(); ++i) {
            const auto nh = nperhashes[i];
            auto &np = nperrows[i];
            if(nh < 3) {
                np = opts.sketchsize_ / nh;
            } else if(nh < 5) {
                np = opts.sketchsize_ / nh * 8;
            } else if(nh < 6) {
                np = opts.sketchsize_ / nh * 6;
            } else np = opts.sketchsize_ / nh * 4;
        }
        auto idx = opts.kmer_result_ > FULL_MMER_SET
        ? SetSketchIndex<uint64_t, LSHIDType>()
        : SetSketchIndex<uint64_t, LSHIDType>(opts.sketchsize_, nperhashes, nperrows);

        // Step 2: Build nearest-neighbor candidate table
        std::vector<std::vector<PairT>> neighbor_lists = build_index(idx, opts, result);
        if(opts.output_kind_ == KNN_GRAPH || opts.output_kind_ == NN_GRAPH_THRESHOLD) {
            refine_results(neighbor_lists, opts, result);
        } else if(opts.output_kind_ == DEDUP) {
            std::fprintf(stderr, "Not yet implemented: de-deduplication via LSH index\n");
        }
        // KNN_GRAPH simply generate top-k
        // 2.
        // NN_GRAPH_THRESHOLD generates all similarities above a threshold
        // 3.
        // DEDUP uses the LSH table to generate only a list of items
    }
}

} // namespace dashing2
