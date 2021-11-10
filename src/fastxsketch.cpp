#include "fastxsketch.h"
#include "mio.hpp"
#include "sketch_core.h"
#include <variant>

//#include <optional>
namespace dashing2 {
using namespace variation;


void FastxSketchingResult::print() {
    std::fprintf(stderr, "%s\n", str().data());
}

using BKRegT = std::conditional_t<(sizeof(RegT) == 4), uint32_t, std::conditional_t<(sizeof(RegT) == 8), uint64_t, u128_t>>;

template<typename C, typename T>
void pop_push(C &c, T &&x, size_t k) {
    if(c.size() < k) c.push(std::move(x));
    else if(x < c.top()) {c.pop(); c.push(std::move(x));}
}

template<typename SrcT, typename CountT=uint32_t>
void bottomk(const std::vector<SrcT> &src, std::vector<BKRegT> &ret, double threshold=0., const CountT *ptr=(CountT *)nullptr, int weighted=-1) {
    if(weighted < 0) weighted = ptr != 0;
    const size_t k = ret.size(), sz = src.size();
    std::priority_queue<BKRegT> pq;
    std::priority_queue<std::pair<double, BKRegT>> wpq;
    for(size_t i = 0; i < sz; ++i) {
        const auto item = src[i];
        const CountT count = ptr ? ptr[i]: CountT(1);
        if(count > threshold) {
            if(weighted) {
                const std::pair<double, BKRegT> key {double(item / count), item};
                pop_push(wpq, key, k);
            } else {
                const BKRegT key = item;
                pop_push(pq, key, k);
            }
        }
    }
    if(weighted) {
        for(size_t i = k; i > 0;ret[--i] = wpq.top().second, wpq.pop());
    } else {
        for(size_t i = k; i > 0;ret[--i] = pq.top(), pq.pop());
    }
}

template<typename T, size_t chunk_size = 65536>
size_t load_copy(const std::string &path, T *ptr, double *cardinality) {
    T *const origptr = ptr;
    if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".gz")) {
        gzFile fp = gzopen(path.data(), "rb");
        if(!fp) return 0; //THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file at ") + path));
        gzread(fp, cardinality, sizeof(*cardinality));
        for(int nr;
            !gzeof(fp) && (nr = gzread(fp, ptr, sizeof(T) * chunk_size)) == sizeof(T) * chunk_size;
            ptr += nr / sizeof(T));
        gzclose(fp);
        return ptr - origptr;
    } else if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".xz")) {
        auto cmd = std::string("xz -dc ") + path;
        std::FILE *fp = ::popen(cmd.data(), "r");
        if(fp == 0) return 0;
        std::fread(cardinality, sizeof(*cardinality), 1, fp);
        for(auto up = (uint8_t *)ptr;!std::feof(fp) && std::fread(up, sizeof(T), chunk_size, fp) == chunk_size; up += chunk_size * sizeof(T));
        ::pclose(fp);
        return ptr - origptr;
    }
    std::FILE *fp = bfopen(path.data(), "rb");
    if(!fp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open ") + path));
    std::fread(cardinality, sizeof(*cardinality), 1, fp);
    const int fd = ::fileno(fp);
    size_t sz = 0;
    if(!::isatty(fd)) {
        struct stat st;
        if(::fstat(fd, &st)) THROW_EXCEPTION(std::runtime_error(std::string("Failed to fstat") + path));
        if(!st.st_size) {
            std::fprintf(stderr, "Warning: Empty file found at %s\n", path.data());
            return 0;
        }
        size_t expected_bytes = st.st_size - 8;
        size_t nb = std::fread(ptr, 1, expected_bytes, fp);
        if(nb != expected_bytes) {
            THROW_EXCEPTION(std::runtime_error("Error in reading from file"));
        }
        sz = expected_bytes / sizeof(T);
    } else {
        auto up = (uint8_t *)ptr;
        for(;!std::feof(fp) && std::fread(up, sizeof(T), chunk_size, fp) == chunk_size; up += chunk_size * sizeof(T));
        sz = (up - (uint8_t *)ptr) / sizeof(T);
    }
    std::fclose(fp);
    return sz;
}

std::string FastxSketchingResult::str() const {
    std::string msg = "FastxSketchingResult @" + to_string(this) + ';';
    if(names_.size()) {
        if(names_.size() < 10) {
            for(const auto &n: names_) msg += n + ",";
        }
        msg += to_string(names_.size()) + " names;";
    }
    if(auto pfsz(nperfile_.size()); pfsz > 0) {
        msg += "sketchedbysequence, ";
        msg += to_string(pfsz) + " seqs";
    } else {msg += "sketchbyline";}
    msg += ';';
    if(signatures_.size()) {
        msg += to_string(signatures_.size()) + " signatures;";
    }
    if(kmers_.size()) {
        msg += to_string(kmers_.size()) + " kmers;";
    }
    if(auto kcsz = kmercounts_.size()) {
        msg += to_string(kcsz) + " kmercounts;";
        long double s = 0., ss = 0.;
        for(const auto v: kmercounts_)
            s += v, ss += v * v;
        msg += "mean: ";
        msg += to_string(double(s / kcsz));
        std::cerr << msg << '\n';
        msg = msg + ", std " + to_string(double(std::sqrt(ss / kcsz - std::pow(s / kcsz, 2.))));
        std::cerr << msg << '\n';
    }
    return msg;
}

INLINE double compute_cardest(const RegT *ptr, const size_t m) {
    double s = 0.;
#if _OPENMP >= 201307L
    #pragma omp simd reduction(+:s)
#endif
    for(size_t i = 0; i < m; ++i) {
        s += ptr[i];
    }
    DBG_ONLY(std::fprintf(stderr, "Sum manually is %g, compared to accumulate with ld %g. diff: %0.20Lg\n", s, double(std::accumulate(ptr, ptr + m, 0.L)), std::accumulate(ptr, ptr + m, 0.L) - static_cast<long double>(s));)
    return m / s;
}




FastxSketchingResult &fastx2sketch(FastxSketchingResult &ret, Dashing2Options &opts, const std::vector<std::string> &paths, std::string outpath) {
    if(paths.empty()) THROW_EXCEPTION(std::invalid_argument("Can't sketch empty path set"));
    std::vector<std::pair<size_t, uint64_t>> filesizes = get_filesizes(paths);
    const size_t nt = std::max(opts.nthreads(), 1u);
    const size_t ss = opts.sketchsize();
    KSeqHolder kseqs(nt);
    std::vector<BagMinHash> bmhs;
    std::vector<ProbMinHash> pmhs;
    std::vector<OPSetSketch> opss;
    std::vector<FullSetSketch> fss;
    std::vector<OrderMinHash> omhs;
    std::vector<Counter> ctrs;
    std::vector<VSetSketch> cfss;
    static_assert(sizeof(pmhs[0].res_[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(bmhs[0].track_ids_[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(opss[0].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(fss[0].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
    auto make = [&](auto &x) {
        x.reserve(nt);
        for(size_t i = 0; i < nt; ++i)
            x.emplace_back(ss);
    };
    auto make_save = [&](auto &x) {
        x.reserve(nt);
        for(size_t i = 0; i < nt; ++i)
            x.emplace_back(ss, opts.save_kmers_, opts.save_kmercounts_);
    };
    if(opts.sspace_ == SPACE_SET) {
        if(opts.kmer_result_ == ONE_PERM) {
            make(opss);
            for(auto &x: opss) x.set_mincount(opts.count_threshold_);
        } else if(opts.kmer_result_ == FULL_SETSKETCH) {
            if(opts.sketch_compressed()) {
                cfss.reserve(nt);
                for(size_t i = 0; i < nt; ++i) {
                    if(opts.fd_level_ == .5) {
                        cfss.emplace_back(NibbleSetS(ss, opts.compressed_b_, opts.compressed_a_));
                    } else if(opts.fd_level_ == 1.) {
                        cfss.emplace_back(ByteSetS(ss, opts.compressed_b_, opts.compressed_a_));
                    } else if(opts.fd_level_ == 2.) {
                        cfss.emplace_back(ShortSetS(ss, opts.compressed_b_, opts.compressed_a_));
                    } else if(opts.fd_level_ == 4.) {
                        cfss.emplace_back(UintSetS(ss, opts.compressed_b_, opts.compressed_a_));
                    }
                }
            } else {
                fss.reserve(nt);
                for(size_t i = 0; i < nt; ++i)
                    fss.emplace_back(opts.count_threshold_, ss, opts.save_kmers_, opts.save_kmercounts_);
            }
        }
    } else if(opts.sspace_ == SPACE_MULTISET) make_save(bmhs);
    else if(opts.sspace_ == SPACE_PSET) make(pmhs);
    else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        if(opts.parse_by_seq_) {
            omhs.reserve(nt);
            for(size_t i = 0; i < nt; omhs.emplace_back(ss, opts.k_), ++i);
        } else {
            THROW_EXCEPTION(std::invalid_argument("Space edit distance is only available in parse-by-seq mode, as it is only defined on strings rather than string collections."));
        }
    }
    while(ctrs.size() < nt) ctrs.emplace_back(opts.cssize());
#define __RESET(tid) do { \
        if(!opss.empty()) opss[tid].reset();\
        else if(!fss.empty()) fss[tid].reset();\
        else if(!cfss.empty()) std::visit([](auto &x) {x.clear();}, cfss[tid]);\
        else if(!bmhs.empty()) bmhs[tid].reset();\
        else if(!pmhs.empty()) pmhs[tid].reset();\
        /*else if(!omhs.empty()) omhs[tid].reset();*/\
        if(ctrs.size() > unsigned(tid)) ctrs[tid].reset();\
    } while(0)

    const uint64_t nitems = paths.size();
    std::string kmeroutpath, kmernamesoutpath;
    if(outpath.size() && outpath != "-" && outpath != "/dev/stdout") {
        const size_t offset = sizeof(nitems) * 2 + sizeof(double) * nitems;
        ::truncate(outpath.data(), offset);
        ret.signatures_.assign(outpath, offset);
        if(opts.save_kmers_) {
            kmeroutpath = outpath + ".kmer64";
            kmernamesoutpath = kmeroutpath + ".names.txt";
        }
    }
    if(kmeroutpath.size()) {
        std::FILE *fp = bfopen(kmeroutpath.data(), "w");
        uint32_t dtype = (uint32_t)opts.input_mode() | (int(opts.canonicalize()) << 8);
        uint32_t sketchsize = opts.sketchsize_;
        uint32_t k = opts.k_;
        uint32_t w = opts.w_ < 0 ? opts.k_: opts.w_;
        checked_fwrite(fp, &dtype, sizeof(dtype));
        checked_fwrite(fp, &sketchsize, sizeof(sketchsize));
        checked_fwrite(fp, &k, sizeof(k));
        checked_fwrite(fp, &w, sizeof(w));
        checked_fwrite(fp, &opts.seedseed_, sizeof(opts.seedseed_));
        if((fp = bfreopen(kmernamesoutpath.data(), "wb", fp)) == 0) THROW_EXCEPTION(std::runtime_error("Failed to open "s + kmernamesoutpath + " for writing."));
        if(bns::filesize(kmeroutpath.data()) != 24) THROW_EXCEPTION(std::runtime_error("kmer out path is the wrong size (expected 16, got "s + std::to_string(bns::filesize(kmeroutpath.data()))));
        static_assert(sizeof(uint32_t) * 4 + sizeof(uint64_t) == 24, "Sanity check");
        ret.kmers_.assign(kmeroutpath, 24);
        for(const auto &n: paths) {
            checked_fwrite(n.data(), 1, n.size(), fp);
            std::fputc('\n', fp);
        }
        std::fclose(fp);
    }
    const int sigshift = opts.sigshift();
    const size_t sigvecsize64 = nitems * ss >> sigshift;
    ret.signatures_.resize(sigvecsize64);
    if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        THROW_EXCEPTION(std::runtime_error("edit distance is only available in parse by seq mode"));
    }
    ret.destination_files_.resize(nitems);
    if(opts.save_kmers_) {
        ret.kmerfiles_.resize(nitems);
    }
    if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
        ret.kmercountfiles_.resize(nitems);
    }
    ret.cardinalities_.resize(nitems, -1.);
#ifndef NDEBUG
    for(size_t i = 0; i < ret.names_.size(); ++i) {
        std::fprintf(stderr, "name %zu is %s\n", i, ret.names_[i].data());
    }
    std::fprintf(stderr, "kmer result type: %s\n", to_string(opts.kmer_result_).data());
    std::fprintf(stderr, "sketching space type: %s\n", to_string(opts.sspace_).data());
#endif
    // We make an exception for iskmer - we only use this if
    //
    if(opts.save_kmers_ && opts.kmer_result_ != FULL_MMER_SEQUENCE) {
        ret.kmers_.resize(ss * nitems);
    }
    if(opts.save_kmercounts_) {
        ret.kmercounts_.resize(ss * nitems);
    }
    if(opts.kmer_result_ == FULL_MMER_SET) {
        ret.kmerfiles_.resize(ret.destination_files_.size());
    }
    OMP_PFOR_DYN
    for(size_t i = 0; i < nitems; ++i) {
        int tid = 0;
        OMP_ONLY(tid = omp_get_thread_num();)
        //const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        //const auto starttime = std::chrono::high_resolution_clock::now();
        auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
        const size_t mss = ss * myind;
        auto &path = paths[myind];
        //std::fprintf(stderr, "parsing from path = %s\n", path.data());
        auto &destination = ret.destination_files_[myind];
        destination = makedest(opts, path, opts.kmer_result_ == FULL_MMER_COUNTDICT);
        const std::string destination_prefix = destination.substr(0, destination.find_last_of('.'));
        std::string kmer_destination_prefix = makedest(opts, path, true);
        kmer_destination_prefix = kmer_destination_prefix.substr(0, kmer_destination_prefix.find_last_of('.'));
        std::string destkmercounts = destination_prefix + ".kmercounts.f64";
        std::string destkmer = kmer_destination_prefix + ".kmer.u64";
        int dkt, dct, dft;
        bool dkif = check_compressed(destkmer, dkt);
        const bool destisfile = check_compressed(destination, dft);
        if(!dkif && opts.kmer_result_ == FULL_MMER_SET && destisfile) {
            dkif = 1; destkmer = destination;
        }
        const bool dkcif = check_compressed(destkmercounts, dct);
        if(ret.kmercountfiles_.size() > myind) ret.kmercountfiles_[myind] = destkmercounts;
        if(opts.cache_sketches_ &&
           (destisfile || (opts.kmer_result_ == FULL_MMER_COUNTDICT && dkif)) &&
           (!opts.save_kmers_ || dkif) &&
           ((!opts.save_kmercounts_ && opts.kmer_result_ != FULL_MMER_COUNTDICT) || dkcif)
        )
        {
            if(opts.kmer_result_ < FULL_MMER_SET) {
                if(ret.signatures_.size()) {
                    if(opts.sketch_compressed()) {
                        std::FILE *ifp = std::fopen(destination.data(), "rb");
                        std::fread(&ret.cardinalities_[myind], sizeof(double), 1, ifp);
                        std::array<long double, 4> arr;
                        std::fread(arr.data(), sizeof(long double), arr.size(), ifp);
                        auto &[a, b, fd_level, sketchsize] = arr;
                        if(fd_level != opts.fd_level_) {
                            std::fprintf(stderr, "fd level found %g mismatches expected %g\n", double(fd_level), opts.fd_level_);
                            THROW_EXCEPTION(std::runtime_error("fd level mismatch."));
                        }
                        if(sketchsize != ss) {
                            std::fprintf(stderr, "fd level found %g mismatches expected %zu\n", double(sketchsize), ss);
                            THROW_EXCEPTION(std::runtime_error("sketch size mismatch."));
                        }
                        RegT *const ptr = &ret.signatures_[(ss >> sigshift) * myind];
                        if(std::fread(ptr, sizeof(RegT), ss >> sigshift, ifp) != (ss >> sigshift)) THROW_EXCEPTION(std::runtime_error("Failed to read compressed signatures from file "s + destination));
                        if(std::fgetc(ifp) != EOF) {
                            THROW_EXCEPTION(std::runtime_error("File corrupted - ifp should be at eof."));
                        }
                        std::fclose(ifp);
                    } else if(load_copy(destination, &ret.signatures_[mss], &ret.cardinalities_[myind]) == 0) {
                        std::fprintf(stderr, "Sketch was not available in file %s... resketching.\n", destination.data());
                        goto perform_sketch;
                    }
                    //ret.cardinalities_[myind] = compute_cardest(&ret.signatures_[mss], ss);
                    DBG_ONLY(std::fprintf(stderr, "Sketch was loaded from %s and has card %g\n", destination.data(), ret.cardinalities_[myind]);)
                }
                if(ret.kmers_.size())
                    load_copy(destkmer, &ret.kmers_[mss], &ret.cardinalities_[myind]);
                if(ret.kmercounts_.size())
                    load_copy(destkmercounts, &ret.kmercounts_[mss], &ret.cardinalities_[myind]);
            } else if(opts.kmer_result_ <= FULL_MMER_SEQUENCE) {
                DBG_ONLY(std::fprintf(stderr, "Cached at path %s, %s, %s\n", destination.data(), destkmercounts.data(), destkmer.data());)
            }
            if(ret.kmerfiles_.size() > myind) {
                ret.kmerfiles_[myind] = destkmer;
            }
            continue;
        } else {
#ifndef NDEBUG
            std::fprintf(stderr, "We skipped caching because with %d as cache sketches\n", opts.cache_sketches_);
            std::fprintf(stderr, "destisfile: %d. is countdict %d. is kmerfile %d\n", destisfile, opts.kmer_result_ == FULL_MMER_COUNTDICT, dkif);
            std::fprintf(stderr, "kc save %d, kmer result %s, dkcif %d\n", opts.save_kmercounts_, to_string(opts.kmer_result_).data(), dkcif);
#endif
        }
        perform_sketch:
        __RESET(tid);
        auto perf_for_substrs = [&](const auto &func) __attribute__((always_inline)) {
            for_each_substr([&](const std::string &subpath) {
                auto lfunc = [&](auto x) __attribute__((always_inline)) {
                    x = maskfn(x);
                    if((!opts.fs_ || !opts.fs_->in_set(x)) && opts.downsample_pass()) func(x);
                };
                auto lfunc2 = [&func](auto x) __attribute__((always_inline)) {func(maskfn(x));};
                const auto seqp = kseqs.kseqs_ + tid;
#define FUNC_FE(f) \
do {\
    if(!opts.fs_ && opts.kmer_downsample_frac_ == 1.) {\
        f(lfunc2, subpath.data(), seqp);\
    } else {\
        f(lfunc, subpath.data(), seqp);\
    } \
} while(0)
                if(opts.use128()) {
                    if(unsigned(opts.k_) <= opts.nremperres128()) {
                        if(entmin) {
                            auto encoder(opts.enc_.to_entmin128());
                            FUNC_FE(encoder.for_each);
                        } else {
                            auto encoder(opts.enc_.to_u128());
                            FUNC_FE(encoder.for_each);
                        }
                    } else {
                        FUNC_FE(opts.rh128_.for_each_hash);
                    }
                } else if(unsigned(opts.k_) <= opts.nremperres64()) {
                    if(entmin) {
                        auto encoder(opts.enc_.to_entmin64());
                        FUNC_FE(encoder.for_each);
                    } else {
                        auto encoder(opts.enc_);
                        FUNC_FE(encoder.for_each);
                    }
                } else {
                    FUNC_FE(opts.rh_.for_each_hash);
                }
#undef FUNC_FE
            }, path);
        };
        if(
            (opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT)
        )
        {
            auto &ctr = ctrs[tid];
            perf_for_substrs([&ctr](auto x) {ctr.add(x);});
            std::vector<u128_t> kmervec128;
            std::vector<uint64_t> kmervec64;
            std::vector<double> kmerveccounts;
            if(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                if(opts.use128())
                    ctr.finalize(kmervec128, kmerveccounts, opts.count_threshold_);
                else
                    ctr.finalize(kmervec64, kmerveccounts, opts.count_threshold_);
                ret.cardinalities_[myind] =
                    opts.kmer_result_ == FULL_MMER_SET ? std::max(kmervec128.size(), kmervec64.size())
                                                       : std::accumulate(kmerveccounts.begin(), kmerveccounts.end(), size_t(0));
            } else if(opts.sspace_ == SPACE_MULTISET) {
                ctr.finalize(bmhs[tid], opts.count_threshold_);
                ret.cardinalities_[myind] = bmhs[tid].total_weight();
                std::copy(bmhs[tid].data(), bmhs[tid].data() + ss, &ret.signatures_[mss]);
            } else if(opts.sspace_ == SPACE_PSET) {
                ctr.finalize(pmhs[tid], opts.count_threshold_);
                std::copy(pmhs[tid].data(), pmhs[tid].data() + ss, &ret.signatures_[mss]);
                ret.cardinalities_[myind] = pmhs[tid].total_weight();
            } else THROW_EXCEPTION(std::runtime_error("Unexpected space for counter-based m-mer encoding"));
                // Make bottom-k if we generated full k-mer sets or k-mer count dictionaries, and copy then over
            if(kmervec64.size() || kmervec128.size()) {
                //std::fprintf(stderr, "If we gathered full k-mers, and we asked for signatures, let's store bottom-k hashes in the signature space\n");
                if(ret.signatures_.size()) {
                    std::vector<BKRegT> keys(ss);
                    double *const kvcp = kmerveccounts.empty() ? static_cast<double *>(nullptr): kmerveccounts.data();
                    if(kmervec128.size()) bottomk(kmervec128, keys, opts.count_threshold_, kvcp);
                    else bottomk(kmervec64, keys, opts.count_threshold_, kvcp);
                    std::copy(keys.begin(), keys.end(), (BKRegT *)&ret.signatures_[mss]);
                }
            }
            std::FILE * ofp = bfopen(destination.data(), "wb");
            checked_fwrite(&ret.cardinalities_[myind], sizeof(ret.cardinalities_[myind]), 1, ofp);
            if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open std::FILE * at") + destination));
            const void *buf = nullptr;
            size_t nb;
            const RegT *srcptr = nullptr;
            if(kmervec128.size()) {
                buf = (const void *)kmervec128.data();
                nb = kmervec128.size() * sizeof(u128_t);
            } else if(kmervec64.size()) {
                buf = (const void *)kmervec64.data();
                nb = kmervec64.size() * sizeof(uint64_t);
            } else if(opts.sspace_ == SPACE_MULTISET) {
                buf = (const void *)bmhs[tid].data();
                nb = ss * sizeof(RegT);
                srcptr = bmhs[tid].data();
            } else if(opts.sspace_ == SPACE_PSET) {
                buf = (const void *)pmhs[tid].data();
                nb = ss * sizeof(RegT);
                srcptr = pmhs[tid].data();
            } else nb = 0, srcptr = nullptr;
            if(srcptr && ret.signatures_.size())
                std::copy(srcptr, srcptr + ss, &ret.signatures_[mss]);
            checked_fwrite(ofp, buf, nb);
            if(opts.save_kmers_ && !(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_COUNTDICT)) {
                assert(ret.kmerfiles_.size());
                ret.kmerfiles_[myind] = destkmer;
                const uint64_t *ptr = opts.sspace_ == SPACE_PSET ? pmhs[tid].ids().data():
                                  opts.sspace_ == SPACE_MULTISET ? bmhs[tid].ids().data():
                                  opts.kmer_result_ == ONE_PERM ? opss[tid].ids().data() :
                                  opts.kmer_result_ == FULL_SETSKETCH ? fss[tid].ids().data():
                                      static_cast<uint64_t *>(nullptr);
                if(!ptr) THROW_EXCEPTION(std::runtime_error("This shouldn't happen"));
                DBG_ONLY(std::fprintf(stderr, "Opening destkmer %s\n", destkmer.data());)
                if((ofp = bfreopen(destkmer.data(), "wb", ofp)) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to write k-mer file"));
                //std::fprintf(stderr, "Writing to file %s\n", destkmer.data());

                checked_fwrite(ofp, ptr, sizeof(uint64_t) * ss);
                DBG_ONLY(std::fprintf(stderr, "About to copy to kmers of size %zu\n", ret.kmers_.size());)
                if(ret.kmers_.size())
                    std::copy(ptr, ptr + ss, &ret.kmers_[mss]);
            }
            if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                //std::fprintf(stderr, "About to save kmer counts manually\n");
                assert(ret.kmercountfiles_.size());
                ret.kmercountfiles_.at(i) = destkmercounts;
                if((ofp = bfreopen(destkmercounts.data(), "wb", ofp)) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to write k-mer counts"));
                std::vector<double> tmp(ss);
#define DO_IF(x) if(x.size()) {std::copy(x[tid].idcounts().begin(), x[tid].idcounts().end(), tmp.data());}
                if(opts.kmer_result_ == FULL_MMER_COUNTDICT || (opts.kmer_result_ == FULL_MMER_SET && opts.save_kmercounts_)) {
                    DBG_ONLY(std::fprintf(stderr, "kvc size %zu. Writing to file %s\n", kmerveccounts.size(), destkmercounts.data());)
                    tmp = std::move(kmerveccounts);
                } else DO_IF(pmhs) else DO_IF(bmhs) else DO_IF(opss) else DO_IF(fss)
#undef DO_IF
                const size_t nb = tmp.size() * sizeof(double);
                checked_fwrite(ofp, tmp.data(), nb);
                if(ret.kmercounts_.size()) {
                    std::fprintf(stderr, "Copying range of size %zu from tmp to ret.kmercounts of size %zu\n", tmp.size(), ret.kmercounts_.size());
                    std::copy(tmp.begin(), tmp.begin() + ss, &ret.kmercounts_[mss]);
                }
            }
            std::fclose(ofp);
        } else if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
            ret.kmers_.clear();
            //std::fprintf(stderr, "Full mmer sequence\n");
            std::FILE * ofp;
            if((ofp = bfopen(destination.data(), "wb")) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open file for writing minimizer sequence"));
            void *dptr = nullptr;
            size_t m = 1 << 18;
            size_t l = 0;
            if(posix_memalign(&dptr, 16, (1 + opts.use128()) * m * sizeof(uint64_t))) THROW_EXCEPTION(std::bad_alloc());

            perf_for_substrs([&](auto x) {
                using DT = decltype(x);
                DT *ptr = (DT *)dptr;
                if(opts.homopolymer_compress_minimizers_ && l > 0 && ptr[l - 1] == x) return;
                if(l == m) {
                    size_t newm = m << 1;
                    void *newptr = nullptr;
                    if(posix_memalign((void **)&newptr, 16, newm * sizeof(DT))) THROW_EXCEPTION(std::bad_alloc());
                    std::copy(ptr, ptr + m, (DT *)newptr);
                    dptr = newptr;ptr = (DT *)dptr;
                    m = newm;
                }
                ptr[l++] = x;
            });
            assert(dptr);
            checked_fwrite(ofp, dptr, l * (1 + opts.use128()) * sizeof(uint64_t));
            ret.cardinalities_[myind] = l;
            std::free(dptr);
            std::fclose(ofp);
        } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
            std::FILE * ofp;
            if((ofp = bfopen(destination.data(), "wb")) == nullptr)
                THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file") + destination + "for writing minimizer sequence"));
            if(opss.empty() && fss.empty() && cfss.empty()) THROW_EXCEPTION(std::runtime_error("Both opss and fss are empty\n"));
            const size_t opsssz = opss.size();
            auto &cret = ret.cardinalities_[myind];
            if(opsssz) {
                //std::fprintf(stderr, "Encode for the opset sketch, %zu is the size for tid %d\n", opss.size(), tid);
                assert(opss.size() > unsigned(tid));
                assert(opss.at(tid).total_updates() == 0);
                auto p = &opss[tid];
                perf_for_substrs([p](auto hv) {p->update(hv);});
                //std::fprintf(stderr, "Encode for the opset sketch. card now: %g, %zu updates\n", opss[tid].getcard(), opss[tid].total_updates());
                assert(ret.cardinalities_.size() > i);
                cret = p->getcard();
            } else {
                //std::fprintf(stderr, "Encode for the set sketch\n");
                if(opts.sketch_compressed()) {
                    std::visit([&](auto &x) {
                        perf_for_substrs([&x](auto hv) {x.update(hv);});
                        cret = x.cardinality();
                        //std::fprintf(stderr, "cret: %g\n", cret);
                    }, cfss[tid]);
                } else {
                    perf_for_substrs([p=&fss[tid]](auto hv) {p->update(hv);});
                    cret = fss[tid].getcard();
                }
            }
            checked_fwrite(ofp, &cret, sizeof(double));
            std::fflush(ofp);
            const uint64_t *ids = nullptr;
            const uint32_t *counts = nullptr;
            // Update this and VSetSketch above to filter down
            const RegT *ptr = opsssz ? opss[tid].data(): fss.size() ? fss[tid].data(): getdata(cfss[tid]);
            const size_t regsize = opsssz  || fss.size() ? sizeof(RegT): size_t(opts.fd_level_);
            assert(ptr);
            if(opts.save_kmers_)
                ids = opsssz ? opss[tid].ids().data(): fss[tid].ids().data();
            if(opts.save_kmercounts_)
                counts = opsssz ? opss[tid].idcounts().data(): fss[tid].idcounts().data();
            if(opts.sketch_compressed()) {
                std::array<long double, 4> arr{opts.compressed_a_, opts.compressed_b_, static_cast<long double>(opts.fd_level_), static_cast<long double>(opts.sketchsize_)};
                checked_fwrite(arr.data(), sizeof(long double), arr.size(), ofp);
                if(opts.fd_level_ == 0.5) {
                    const uint8_t *srcptr = std::get<NibbleSetS>(cfss[tid]).data();
                    for(size_t i = 0; i < opts.sketchsize_; i += 2) {
                        uint8_t reg = (srcptr[i] << 4) | srcptr[i + 1];
                        checked_fwrite(ptr, sizeof(reg), 1, ofp);
                    }
                } else {
                    checked_fwrite(ptr, sizeof(RegT), ss >> sigshift, ofp);
                }
            } else {
                ::write(::fileno(ofp), ptr, regsize * ss);
            }
            std::fclose(ofp);
            if(ptr && ret.signatures_.size()) {
                if(opts.fd_level_ != .5) {
                    std::copy(ptr, ptr + (ss >> sigshift), &ret.signatures_[mss >> sigshift]);
                } else {
                    const uint8_t *srcptr = std::get<NibbleSetS>(cfss[tid]).data();
                    uint8_t *destptr = (uint8_t *)&ret.signatures_[mss >> sigshift];
                    for(size_t i = 0; i < opts.sketchsize_; i += 2) {
                        *destptr++ = (srcptr[i] << 4) | srcptr[i + 1];
                    }
                }
            }
            if(ids && ret.kmers_.size())
                std::copy(ids, ids + ss, &ret.kmers_[mss]);
            if(counts && ret.kmercounts_.size())
                std::copy(counts, counts + ss, &ret.kmercounts_[mss]);
        } else THROW_EXCEPTION(std::runtime_error("Unexpected: Not FULL_MMER_SEQUENCE, FULL_MMER_SET, ONE_PERM, FULL_SETSKETCH, SPACE_MULTISET, or SPACE_PSET"));
    } // parallel paths loop
    ret.names_ = paths;
    return ret;
}




} // dashing2
