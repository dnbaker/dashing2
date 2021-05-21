#include "cmp_main.h"
#include "sketch/count_eq.h"
namespace dashing2 {
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
                std::fprintf(stderr, "Double reg = %u\n", ptr[i], (reg2sig(sig1) & 0xfu), (reg2sig(sig2) & 0xfu) << 4);
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
        std::fprintf(stderr, "fd: %g\n", fd);
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
            uint8_t maxreg = 0, minreg = 15;
            uint8_t *ptr = (uint8_t *)compressed_reps;
            for(size_t i = 0; i < nsigs / 2; ++i) {
                uint8_t lower_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i] / a) * logbinv))));
                uint8_t upper_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i + 1] / a) * logbinv)))) << 4;
                maxreg = std::max(lower_half, std::max(uint8_t(upper_half >> 4), maxreg));
                minreg = std::min(lower_half, std::min(uint8_t(upper_half >> 4), minreg));
                ptr[i] = lower_half | upper_half;
                std::fprintf(stderr, "lower half: %d. upper half: %d. combined: %d\n", lower_half, upper_half, ptr[i]);
            }
            std::fprintf(stderr, "Smallest register %d, max %d\n", maxreg, minreg);
        }
    }
    return ret;
}
void emit_all_pairs(Dashing2DistOptions &opts, const SketchingResult &result) {
    const size_t ns = result.names_.size();
#ifndef NDEBUG
    for(size_t i = 0; i < ns; ++i) {
        std::fprintf(stderr, "name %s for index %zu\n", result.names_[i].data(), i);
    }
#endif
    auto compare = [&](size_t i, size_t j) -> double {
        double ret = std::numeric_limits<double>::max();
        const double b2pow = -std::ldexp(1., -static_cast<int>(opts.fd_level_ * 8.));
        const double ib2pow = 1. / (1. + b2pow);
        const double invdenom = 1. / opts.sketchsize_;
        const double lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
        const double poisson_mult = 1. / std::max(1, opts.k_);
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
                        std::fprintf(stderr, "gt/lt/eq: %zu/%zu/%zu\n", res.first, res.second, opts.sketchsize_ - res.first - res.second);
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
                const double b = opts.compressed_b_;
                double mu;
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
                double alpha = gtlt.first * invdenom;
                double beta = gtlt.second * invdenom;
                double eq = (1. - alpha - beta);
                double lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
                double ucard = alpha + beta >= 1. ? lhcard + rhcard: std::max((lhcard + rhcard) / (2. - alpha - beta), 0.);
                std::fprintf(stderr, "alpha %g, beta %g, eq = %g, lhcard %g, rhcard %g, ucard %g\n", alpha, beta, eq, lhcard, rhcard, ucard);
                std::fprintf(stderr, "gtlt %zu/%zu\n", gtlt.first, gtlt.second);
                double isz = ucard * eq;
                double sim = isz / ucard;
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
    };
    std::FILE *ofp = opts.outfile_path_.empty() ? stdout: std::fopen(opts.outfile_path_.data(), "wb");
    if(!ofp) throw std::runtime_error(std::string("Failed to open path at ") + opts.outfile_path_);
    std::setvbuf(ofp, nullptr, _IOFBF, 1<<17);
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::thread sub;
    for(size_t i = 0; i < ns; ++i) {
        size_t nelem = asym ? ns: ns - i - 1;
        std::unique_ptr<float[]> dat(new float[nelem]);
        const size_t diff = asym ? 0: i + 1;
        const size_t nwritten = asym ? ns: ns - i - 1;
        OMP_PFOR_DYN
        for(size_t start = asym ? 0: i + 1;start < ns; ++start) {
            auto v = compare(i, start);
            std::fprintf(stderr, "Calling %zu, %zu (i < ns = %zu)\n", i, start, ns);
            dat[start - diff] = v;
        }
        if(sub.joinable()) sub.join();
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
        // Build LSH index, then generate values
        // 1.
        // KNN_GRAPH simply generate top-k
        // 2.
        // NN_GRAPH_THRESHOLD generates all similarities above a threshold
        // 3.
        // DEDUP uses the LSH table to generate only a list of items
    }
}

} // namespace dashing2
