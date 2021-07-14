#include "cmp_main.h"
#include "sketch/count_eq.h"
#include "index_build.h"
#include "refine.h"
#include "emitnn.h"
#include "mio.hpp"
#include "wcompare.h"
#include "levenshtein-sse.hpp"


namespace dashing2 {
//using sketch::lsh::SetSketchIndex;
static INLINE uint64_t reg2sig(RegT x) {
    uint64_t seed = 0;
    CONST_IF(sizeof(RegT) <= 8) {
        std::memcpy(&seed, &x, sizeof(x));
        return wy::wyhash64_stateless(&seed);
    } else {
        std::memcpy(&seed, &x, std::min(sizeof(x), sizeof(seed)));
        uint64_t nextseed = wy::wyhash64_stateless(&seed);
        nextseed ^= ((uint64_t *)&x)[1];
        return wy::wyhash64_stateless(&nextseed);
    }
}

#ifdef _OPENMP
#define OMP_STATIC_SCHED32 _Pragma("omp parallel for schedule(static, 32)")
#else
#define OMP_STATIC_SCHED32
#endif

std::string path2cmd(const std::string &path) {
    if(std::equal(path.rbegin(), path.rbegin() + 3, "zg.")) return std::string("gzip -dc ") + path;
    if(std::equal(path.rbegin(), path.rbegin() + 3, "zx.")) return std::string("xz -dc ") + path;
    if(std::equal(path.rbegin(), path.rbegin() + 4, "2zb.")) return std::string("bzip2 -dc ") + path;
    return std::string("cat ") + path;
}

std::tuple<void *, double, double>
make_compressed(int truncation_method, double fd, const std::vector<RegT> &sigs, bool is_edit_distance) {
    std::tuple<void *, double, double> ret = {nullptr, 0., 0.};
    if(fd >= sizeof(RegT)) return ret;
    size_t mem = fd * sigs.size();
    if(fd == 0. && std::fmod(fd * sigs.size(), 1.)) THROW_EXCEPTION(std::runtime_error("Can't do nibble registers without an even number of signatures"));
    auto &compressed_reps = std::get<0>(ret);
    if(posix_memalign((void **)&compressed_reps, 64, mem)) THROW_EXCEPTION(std::bad_alloc());
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
    } else if(truncation_method > 0) {
        std::fprintf(stderr, "Performing %d-bit compression\n", int(fd * 8.));
        if(fd == 0.5) {
            uint8_t *ptr = (uint8_t *)compressed_reps;
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs / 2; ++i) {
                auto sig1 = sigs[2 * i], sig2 = sigs[2 * i + 1];
                ptr[i] = (reg2sig(sig1) & 0xfu) | ((reg2sig(sig2) & 0xfu) << 4);
            }
        } else {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs; ++i) {
                const auto sig = reg2sig(sigs[i]);
                if(fd == 8) ((uint64_t *)compressed_reps)[i] = sig;
                else if(fd == 4) ((uint32_t *)compressed_reps)[i] = sig;
                else if(fd == 2) ((uint16_t *)compressed_reps)[i] = sig;
                else ((uint8_t *)compressed_reps)[i] = sig;
            }
        }
    } else {
        RegT minreg = sigs[0], maxreg = minreg;
        OMP_ONLY(_Pragma("omp parallel for simd reduction(min:minreg) reduction(max:maxreg)"))
        for(size_t i = 0; i < nsigs; ++i) {
            if(sigs[i] == std::numeric_limits<RegT>::max())
                continue;
            minreg = std::min(minreg, sigs[i]);
            maxreg = std::max(maxreg, sigs[i]);
        }
        double q = fd == 1. ? 254.3: fd == 2. ? 65534.3: fd == 4. ? 4294967294.3: fd == 8. ? 18446744073709551615. : fd == 0.5 ? 15.4: -1.;
        assert(q > 0.);
        auto [b, a] = sketch::CSetSketch<RegT>::optimal_parameters(minreg, maxreg, q);
        std::fprintf(stderr, "Truncated via setsketch, a = %Lg and b = %Lg\n", a, b);
        std::get<1>(ret) = a;
        std::get<2>(ret) = b;
        const RegT logbinv = 1. / std::log1p(b - 1.);
        if(fd == 0.5) {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs / 2; ++i) {
                const uint8_t lower_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i] / a) * logbinv))));
                const uint8_t upper_half = std::max(0, std::min(int(q) + 1, static_cast<int>((1. - std::log(sigs[2 * i + 1] / a) * logbinv)))) << 4;
                ((uint8_t *)compressed_reps)[i] = lower_half | upper_half;
            }
        } else {
            OMP_STATIC_SCHED32
            for(size_t i = 0; i < nsigs; ++i) {
                RegT x = sigs[i];
                const RegT sub = 1. - std::log(x / a) * logbinv;
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
    }
    return ret;
}

LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j) {
    LSHDistType ret = std::numeric_limits<LSHDistType>::max();
    const LSHDistType lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
    const LSHDistType invdenom = 1. / opts.sketchsize_;
    auto sim2dist = [poisson_mult=-1. / std::max(1, opts.k_)](auto x) -> double {if(x) return std::log(2. * x / (1. + x)) * poisson_mult; return std::numeric_limits<double>::infinity();};
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
case v: {\
    VERBOSE_ONLY(std::fprintf(stderr, "Doing comparing between %zu and %zu with %d bits\n", i, j, int(8. * opts.fd_level_));)\
    TYPE *ptr = static_cast<TYPE *>(opts.compressed_ptr_);\
    res = sketch::eq::count_gtlt(ptr + i * opts.sketchsize_, ptr + j * opts.sketchsize_, opts.sketchsize_);} break;
                CASE_ENTRY(16, uint64_t)
                CASE_ENTRY(8, uint32_t)
                CASE_ENTRY(4, uint16_t)
                CASE_ENTRY(2, uint8_t)
#undef CASE_ENTRY
                case 1: {
                    uint8_t *ptr = static_cast<uint8_t *>(opts.compressed_ptr_);
                    res = sketch::eq::count_gtlt_nibbles(ptr + i * opts.sketchsize_ / 2, ptr + j * opts.sketchsize_ / 2, opts.sketchsize_);
                    //std::fprintf(stderr, "gt/lt/eq: %zu/%zu/%zu\n", size_t(res.first), size_t(res.second), size_t(opts.sketchsize_ - res.first - res.second));
                    break;
                }
                default: __builtin_unreachable();
            }
        }
        if(bbit_c) {
            // ret = ((num / denom) - (1. / 2^b)) / (1. - 1. / 2^b);
            // maps equality to 1 and down-estimates for account for collisions
            const LSHDistType b2pow = -std::ldexp(1., -static_cast<int>(opts.fd_level_ * 8.));
            ret = std::max(0., std::fma(res.first, invdenom, b2pow) / (1. + b2pow));
            if(opts.measure_ == INTERSECTION)
                ret *= std::max((lhcard + rhcard) / (2. - (1. - ret)), 0.);
            else if(opts.measure_ == CONTAINMENT)
                ret = std::max((lhcard + rhcard) / (2. - (1. - ret)), 0.) * ret / lhcard;
            else if(opts.measure_ == POISSON_LLR)
                ret = sim2dist(ret);
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
                ret = sim2dist(ret);
        }
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE && (opts.exact_kmer_dist_ || opts.measure_ == M_EDIT_DISTANCE)) {
        //std::fprintf(stderr, "Return exact distance for edit distance between sequences\n");
        assert(result.sequences_.size() > std::max(i, j) || !std::fprintf(stderr, "Expected sequences to be non-null for exact edit distance calculation (%zu vs %zu/%zu)\n", result.sequences_.size(), i, j));
        return levenshteinSSE::levenshtein(result.sequences_[i], result.sequences_[j]);
    } else if(opts.kmer_result_ <= FULL_SETSKETCH) {
        const RegT *lhsrc = &result.signatures_[opts.sketchsize_ * i], *rhsrc = &result.signatures_[opts.sketchsize_ * j];
        if(opts.sspace_ == SPACE_SET) {
            auto gtlt = sketch::eq::count_gtlt(lhsrc, rhsrc, opts.sketchsize_);
            LSHDistType alpha, beta, eq, lhcard, ucard, rhcard;
            alpha = gtlt.first * invdenom;
            beta = gtlt.second * invdenom;
            lhcard = result.cardinalities_.at(i), rhcard = result.cardinalities_.at(j);
            eq = (1. - alpha - beta);
            if(eq <= 0.)
                return opts.measure_ != POISSON_LLR ? 0.: std::numeric_limits<double>::max();
            ucard = std::max(lhcard + rhcard / (2. - alpha - beta), 0.);
            LSHDistType isz = ucard * eq, sim = eq;
            ret = opts.measure_ == SIMILARITY ? sim
                : opts.measure_ == INTERSECTION ? isz
                : opts.measure_ == SYMMETRIC_CONTAINMENT ? isz / (std::min(lhcard, rhcard))
                : opts.measure_ == POISSON_LLR ? sim2dist(sim): LSHDistType(-1);
            assert(ret >= 0. || !std::fprintf(stderr, "measure: %s. sim: %g. isz: %g\n", to_string(opts.measure_).data(), sim, isz));
        } else {
            const auto neq = sketch::eq::count_eq(&result.signatures_[opts.sketchsize_ * i], &result.signatures_[opts.sketchsize_ * j], opts.sketchsize_);
            ret = invdenom * neq;
            DBG_ONLY(std::fprintf(stderr, "Computing number of equal registers between %zu and %zu, resulting in %zu/%zu (%g)\n", i, j, size_t(neq), opts.sketchsize_, ret);)
            if(opts.measure_ == INTERSECTION) {
                ret *= std::max((lhcard + rhcard) / (1. + ret), 0.);
            } else if(opts.measure_ == SYMMETRIC_CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1. + ret), 0.) / std::min(lhcard, rhcard);
            else if(opts.measure_ == CONTAINMENT) ret *= std::max((lhcard + rhcard) / (1. + ret), 0.) / lhcard;
            else if(opts.measure_ == POISSON_LLR) ret = sim2dist(ret);
        }
    } else {
#define CORRECT_RES(res, measure, lhc, rhc)\
            if(measure == INTERSECTION) {\
                /* do nothing*/ \
            } else if(measure == SYMMETRIC_CONTAINMENT) {\
                res = res / std::min(lhc, rhc);\
            } else if(measure == POISSON_LLR || measure == SIMILARITY) {\
                res = res / (lhc + rhc - res);\
            } else if(measure == CONTAINMENT) {\
                res /= lhc;\
            }\
            ret = res;
        const std::string &lpath = result.destination_files_[i], &rpath = result.destination_files_[j];
        if(lpath.empty() || rpath.empty()) THROW_EXCEPTION(std::runtime_error("Destination files for k-mers empty -- cannot load from disk"));
        //std::fprintf(stderr, "Paths %zu/%zu are %s/%s\n", i, j, lpath.data(), rpath.data());
        std::string dummy;
        auto &lnpath(result.kmercountfiles_.size() ? result.kmercountfiles_[i]: dummy);
        auto &rnpath(result.kmercountfiles_.size() ? result.kmercountfiles_[j]: dummy);
        const bool lcomp = iscomp(lpath), rcomp = iscomp(rpath);
        const bool lncomp = iscomp(lnpath), rncomp = iscomp(rnpath);
        if(lcomp || rcomp || lncomp || rncomp)
        {
            std::FILE *lhk = 0, *rhk = 0, *lhn = 0, *rhn = 0;
            std::string lcmd = path2cmd(lpath);
            std::string rcmd = path2cmd(rpath);
            std::cerr << lcmd << '\n';
            std::cerr << rcmd << '\n';
            if((lhk = ::popen(lcmd.data(), "r")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to run lcmd '") + lcmd + "'"));
            if((rhk = ::popen(rcmd.data(), "r")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to run rcmd '") + rcmd + "'"));
            if(result.kmercountfiles_.size()) {
                lcmd = path2cmd(result.kmercountfiles_[i]);
                rcmd = path2cmd(result.kmercountfiles_[j]);
                if((lhn = ::popen(lcmd.data(), "r")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to run lcmd '") + lcmd + "'"));
                if((rhn = ::popen(rcmd.data(), "r")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to run lcmd '") + rcmd + "'"));
            }
            const auto lhc = result.cardinalities_[i], rhc = result.cardinalities_[j];
            if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                if(opts.exact_kmer_dist_) {
                    auto [edit_dist, max_edit_dist] = mmer_edit_distance(lhk, rhk, opts.use128());
                    ret = opts.measure_ == M_EDIT_DISTANCE ? edit_dist: max_edit_dist - edit_dist;
                } else {
                    ret = hamming_compare_f64(lhk, rhk);
                }
            } else {
                std::pair<double, double> wcret = weighted_compare(lhk, rhk, 0, 0, lhc, rhc);
                auto [isz_size, union_size] = wcret;
                double res = isz_size;
                CORRECT_RES(res, opts.measure_, lhc, rhc)
            }
            ::pclose(lhk); ::pclose(rhk);
            if(lhn) ::pclose(lhn), ::pclose(rhn);
        } else {
            //std::fprintf(stderr, "Using mmapping, %s and %s\n", lpath.data(), rpath.data());
            mio::mmap_source lhs(lpath), rhs(rpath);
            std::unique_ptr<mio::mmap_source> lhn, rhn;
            if(result.kmercountfiles_.size()) {
                lhn.reset(new mio::mmap_source(result.kmercountfiles_[i]));
                rhn.reset(new mio::mmap_source(result.kmercountfiles_[j]));
            }
            if(opts.use128()) THROW_EXCEPTION(std::runtime_error("Not yet implemented: 128-bit exact k-mer comparisons"));
            const uint64_t *lptr = (const uint64_t *)lhs.data(), *rptr = (const uint64_t *)rhs.data();
            const size_t lhl = lhs.size() / 8, rhl = rhs.size() / 8;
            if(lhn && rhn) {
                const double *lnptr = (const double *)lhn->data(), *rnptr = (const double *)rhn->data();
                double lhc = result.cardinalities_[i], rhc = result.cardinalities_[j];
                auto [isz_size, union_size] = weighted_compare(lptr, lnptr, lhl, lhc, rptr, rnptr, rhl, rhc);
                double res = isz_size;
                CORRECT_RES(res, opts.measure_, lhc, rhc)
            } else {
                double res = set_compare(lptr, lhl, rptr, rhl);
                CORRECT_RES(res, opts.measure_, lhl, rhl)
            }
        }
#undef CORRECT_RES
        // Compare exact representations, not compressed shrunk
    }
    if(std::isnan(ret) || std::isinf(ret)) ret = std::numeric_limits<double>::max();
    return ret;
}
void emit_all_pairs(Dashing2DistOptions &opts, const SketchingResult &result) {
    std::fprintf(stderr, "Beginning emit_all_pairs\n");
    const size_t ns = result.names_.size();
    std::FILE *ofp = opts.outfile_path_.empty() || opts.outfile_path_ == "-" ? stdout: std::fopen(opts.outfile_path_.data(), "w");
    if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open path at ") + opts.outfile_path_));
    std::setvbuf(ofp, nullptr, _IOFBF, 1<<17);
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::deque<std::pair<std::unique_ptr<float[]>, size_t>> datq;
    volatile int loopint = 0;
    std::mutex datq_lock;
    // Emit Header
    if(opts.output_format_ == HUMAN_READABLE) {
        std::fprintf(ofp, "#Dashing2 %s Output\n", asym ? "Asymmetric pairwise": "PHYLIP pairwise");
        std::fprintf(ofp, "#Dashing2Options: %s\n", opts.to_string().data());
        std::fprintf(ofp, "#Sources");
        for(size_t i = 0; i < result.names_.size(); ++i) {
            std::fprintf(ofp, "\t%s", result.names_[i].data());
        }
        std::fputc('\n', ofp);
        std::fprintf(ofp, "%zu\n", ns);
    }
    std::thread sub = std::thread([&](){
        while(loopint == 0 || datq.size()) {
            if(datq.empty()) {
                std::this_thread::sleep_for(std::chrono::duration<size_t, std::nano>(2000));
            } else {
                auto &f = datq.front();
                const size_t nwritten = asym ? ns: ns - f.second - 1;
                if(opts.output_format_ == HUMAN_READABLE) {
                    auto datp = f.first.get();
                    auto &seqn = result.names_[f.second];
                    std::string fn(seqn);
                    if(fn.size() < 9) fn.append(9 - fn.size(), ' ');
                    std::fwrite(fn.data(), 1, fn.size(), ofp);
                    for(size_t i = 0; i < nwritten;std::fprintf(ofp, "\t%0.9g", datp[i++]));
                    std::fputc('\n', ofp);
                } else if(opts.output_format_ == MACHINE_READABLE) {
                    if(std::fwrite(datq.front().first.get(), sizeof(float), nwritten, ofp) != nwritten)
                        THROW_EXCEPTION(std::runtime_error(std::string("Failed to write row ") + std::to_string(datq.front().second) + " to disk"));
                } // else throw std::runtime_error("This should never happen");
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.pop_front();
            }
        }
    });
    for(size_t i = 0; i < ns; ++i) {
        // TODO: batch queries together for cache efficiency (distmat::parallel_fill for an example)
        size_t nelem = asym ? ns: ns - i - 1;
        std::unique_ptr<float[]> dat(new float[nelem]);
        const auto datp = dat.get() - (asym ? size_t(0): i + 1);
        OMP_PFOR_DYN
        for(size_t start = asym ? 0: i + 1;start < ns; ++start) {
            datp[start] = compare(opts, result, i, start);
        }
        std::lock_guard<std::mutex> guard(datq_lock);
        datq.emplace_back(std::pair<std::unique_ptr<float[]>, size_t>{std::move(dat), i});
    }
    loopint = 1;
    if(sub.joinable()) sub.join();
    if(ofp != stdout) std::fclose(ofp);
}
void cmp_core(Dashing2DistOptions &opts, SketchingResult &result) {
    if(opts.kmer_result_ >= FULL_MMER_SET && result.cardinalities_.size() && result.cardinalities_.front() < 0.) {
        std::fprintf(stderr, "Getting cardinalities from k-mer set/seq/countdict\n");
        const size_t n = result.cardinalities_.size();
        if(opts.parse_by_seq_)
            THROW_EXCEPTION(std::runtime_error("Not yet supported"));
        OMP_PFOR
        for(size_t i = 0; i < n; ++i) {
            int ft;
#if 0
            if(!check_compressed(result.names_.at(i), ft)) {
                std::fprintf(stderr, "Missing input file %s at %zu, ignoring...\n", &result.names_.at(i)[0], i);
            }
#endif
            std::FILE *ifp = nullptr;
            std::string fn = opts.kmer_result_ == FULL_MMER_SET ? result.kmerfiles_.at(i): result.destination_files_.at(i);
            if(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                if(!check_compressed(fn, ft)) throw std::runtime_error(std::string("Missing kmerfile or destination file: ") + fn);
                if(endswith(fn, ".xz")) ft = 2;
                else if(endswith(fn, ".gz")) ft = 1;
                else ft = 0;
                std::string cmd = std::string(ft == 0 ? "cat ": ft == 1 ? "gzip -dc ": ft == 2 ? "xz -dc " : "unknowncommand") + fn;
                std::fprintf(stderr, "Checking with command = '%s'\n", cmd.data());
                if(!(ifp = ::popen(cmd.data(), "r")))
                    THROW_EXCEPTION(std::runtime_error(std::string("Command ") + "'" + cmd + "' failed."));
                size_t c = 0;
                for(uint64_t x;std::fread(&x, sizeof(x), 1, ifp) == 1u;++c)
                result.cardinalities_[i] = c;
            } else if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                if(!check_compressed(result.kmercountfiles_[i], ft)) throw std::runtime_error("Missing kmercountfile");
                if(endswith(result.kmercountfiles_[i], ".xz")) ft = 2;
                else if(endswith(result.kmercountfiles_[i], ".gz")) ft = 1;
                else ft = 0;
                std::string cmd = std::string(ft == 0 ? "cat ": ft == 1 ? "gzip -dc ": ft == 2 ? "xz -dc " : "unknowncommand");
                if(cmd == "unknowncommand") THROW_EXCEPTION(std::runtime_error("Failure"));
                cmd += result.kmercountfiles_[i];
                std::fprintf(stderr, "Calling %s\n", cmd.data());
                if(!(ifp = ::popen(cmd.data(), "r")))
                    THROW_EXCEPTION(std::runtime_error(std::string("Command ") + "'" + cmd + "' failed."));
                double x, c, s;
                for(x = c = s = 0.;std::fread(&x, sizeof(x), 1, ifp) == 1u;sketch::kahan::update(s, c, x));
                result.cardinalities_[i] = s;
            }
            if(ifp) ::pclose(ifp);
        }
    }
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
    std::tie(opts.compressed_ptr_, opts.compressed_a_, opts.compressed_b_) = make_compressed(opts.truncation_method_, opts.fd_level_, result.signatures_, opts.sspace_ == SPACE_EDIT_DISTANCE);
    if(opts.output_kind_ <= ASYMMETRIC_ALL_PAIRS) {
        emit_all_pairs(opts, result);
        return;
    }
    // This is LSH-index assisted KNN graphs +
    // thresholded nn graphs

    // Step 1: Build LSH Index
    std::vector<uint64_t> nperhashes{1, 2, 3, 4, 6, 8};
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
    auto idx = [&]() -> SetSketchIndex<uint64_t, LSHIDType> {
        if(opts.kmer_result_ >= FULL_MMER_SET) return {};
        return SetSketchIndex<uint64_t, LSHIDType>(opts.sketchsize_, nperhashes, nperrows);
    }();

    // Step 2: Build nearest-neighbor candidate table
    if(opts.output_kind_ == KNN_GRAPH || opts.output_kind_ == NN_GRAPH_THRESHOLD) {
        std::vector<pqueue> neighbor_lists = build_index(idx, opts, result);
        refine_results(neighbor_lists, opts, result);
        emit_neighbors(neighbor_lists, opts, result);
    } else if(opts.output_kind_ == DEDUP) {
        std::vector<size_t> ids;
        std::vector<std::vector<size_t>> constituents;
        //const double simthres = 0.9; // This will be changed in the future, but for now, this is hard-coded and dumb
        std::vector<size_t> order(result.names_.size());
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&result](auto x, auto y) {return result.cardinalities_[x] < result.cardinalities_[y];});
        // General strategy:
        // Use a given similarity threshold to then group items into the cluster
        // to which they are most similar if they are > than
        for(size_t idx = 0; idx < order.size(); ++idx) {
            //auto oid = order[idx];
        }
        THROW_EXCEPTION(std::runtime_error("Not implemented: Deduplication"));
    }
}

} // namespace dashing2
