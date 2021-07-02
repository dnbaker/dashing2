#include "fastxsketch.h"
#include "mio.hpp"
#include "sketch_core.h"
//#include <optional>
namespace dashing2 {


#define checked_fwrite(fp, ptr, nb) \
    do {\
        if(unsigned long long lrc = std::fwrite(static_cast<const void *>(ptr), 1, nb, fp); lrc != static_cast<size_t>(nb)) \
            throw std::runtime_error(std::string("[E:") + __PRETTY_FUNCTION__ + ':' + __FILE__ + std::to_string(__LINE__) + "] Failed to perform buffered write of " + std::to_string(static_cast<size_t>(nb)) + " bytes, instead writing " + std::to_string(lrc) + " bytes");\
    } while(0)

void FastxSketchingResult::print() {
    std::fprintf(stderr, "%s\n", str().data());
}

using BKRegT = std::conditional_t<(sizeof(RegT) == 4), uint32_t, std::conditional_t<(sizeof(RegT) == 8), uint64_t, u128_t>>;

template<typename SrcT, typename CountT=uint32_t>
void bottomk(const std::vector<SrcT> &src, std::vector<BKRegT> &ret, double threshold=0., const CountT *ptr=(CountT *)nullptr) {
    const size_t k = ret.size(), sz = src.size();
    std::priority_queue<BKRegT> pq;
    for(size_t i = 0; i < sz; ++i) {
        const auto item = src[i];
        const CountT count = ptr ? ptr[i]: CountT(1);
        if(count > threshold) {
            const BKRegT key = item;
            if(pq.size() < k) pq.push(key);
            else if(key < pq.top()) {
                pq.pop(); pq.push(key);
            }
        }
    }
    for(size_t i = k; i > 0;ret[--i] = pq.top(), pq.pop());
}

template<typename T, size_t chunk_size = 65536>
void load_copy(const std::string &path, T *ptr) {
    if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".gz")) {
        gzFile fp = gzopen(path.data(), "rb");
        if(!fp) throw std::runtime_error(std::string("Failed to open file at ") + path);
        for(int nr;
            !gzeof(fp) && (nr = gzread(fp, ptr, sizeof(T) * chunk_size)) == sizeof(T) * chunk_size;
            ptr += nr / sizeof(T));
        gzclose(fp);
        return;
    } else if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".xz")) {
        std::FILE *fp = ::popen((std::string("xz -dc ") + path).data(), "r");
        for(auto up = (uint8_t *)ptr;!std::feof(fp) && std::fread(up, sizeof(T), chunk_size, fp) == chunk_size; up += chunk_size * sizeof(T));
        ::pclose(fp);
    }
    std::FILE *fp = std::fopen(path.data(), "rb");
    if(!fp) throw std::runtime_error(std::string("Failed to open ") + path);
    struct stat st;
    if(::fstat(::fileno(fp), &st))
        throw std::runtime_error(std::string("Failed to get fd from ") + path + std::to_string(::fileno(fp)));
    DBG_ONLY(std::fprintf(stderr, "Size of file at %s: %zu\n", path.data(), size_t(st.st_size));)
    size_t nb = std::fread(ptr, 1, st.st_size, fp);
    if(nb != size_t(st.st_size)) {
        std::fprintf(stderr, "Failed to copy to ptr %p. Expected to read %zu, got %zu\n", (void *)ptr, size_t(st.st_size), nb);
        throw std::runtime_error("Error in reading from file");
    }
    std::fclose(fp);
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

struct KSeqHolder {
    kseq_t *kseqs_;
    size_t n_;
    KSeqHolder(size_t n): kseqs_(static_cast<kseq_t *>(std::calloc(n, sizeof(kseq_t)))), n_(n) {
        if(!kseqs_) throw std::bad_alloc();
        for(auto p = kseqs_; p < kseqs_ + n_; ++p)
            ks_resize(&p->seq, 1<<20);
    }
    void free_item(kseq_t &seq) {
        std::free(seq.name.s);
        std::free(seq.comment.s);
        std::free(seq.seq.s);
        std::free(seq.qual.s);
        ks_destroy(seq.f);
    }
    ~KSeqHolder() {
        for(size_t i = 0; i < n_; free_item(kseqs_[i++]));
        std::free(kseqs_);
    }
};
INLINE double compute_cardest(const RegT *ptr, const size_t m) {
    double s = 0.;
#if _OPENMP >= 201307L
    #pragma omp simd reduction(+:s)
#endif
    for(size_t i = 0; i < m; ++i)
        s += ptr[i];
    return m / s;
}
static inline bool check_compressed(std::string &path) {
    if(bns::isfile(path))
        return true;
    if(bns::isfile(path + ".gz")) {
        path = path + ".gz";
        return true;
    }
    if(bns::isfile(path + ".xz")) {
        path = path + ".xz";
        return true;
    }
    return false;
}

FastxSketchingResult fastx2sketch(Dashing2Options &opts, const std::vector<std::string> &paths) {
    if(paths.empty()) throw std::invalid_argument("Can't sketch empty path set");
    std::vector<std::pair<size_t, uint64_t>> filesizes = get_filesizes(paths);
    const size_t nt = std::max(opts.nthreads(), 1u);
    const size_t ss = opts.sketchsize();
    FastxSketchingResult ret;
    ret.options_ = &opts;
    std::vector<BagMinHash> bmhs;
    std::vector<ProbMinHash> pmhs;
    std::vector<OPSetSketch> opss;
    std::vector<FullSetSketch> fss;
    std::vector<OrderMinHash> omhs;
    std::vector<Counter> ctrs;
    static_assert(sizeof(pmhs[0].res_[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(bmhs[0].track_ids_[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(opss[0].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
    static_assert(sizeof(fss[0].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
    KSeqHolder kseqs(nt);
    auto make = [&](auto &x) {
        x.reserve(nt);
        for(size_t i = 0; i < nt; ++i)
            x.emplace_back(ss);
    };
    auto make_save = [&](auto &x) {
        x.reserve(nt);
        for(size_t i = 0; i < nt; ++i)
            x.emplace_back(ss, opts.save_kmers_ || opts.build_mmer_matrix_, opts.save_kmercounts_ || opts.build_count_matrix_);
    };
    if(opts.sspace_ == SPACE_SET) {
        if(opts.kmer_result_ == ONE_PERM) {
            make(opss);
            for(auto &x: opss) x.set_mincount(opts.count_threshold_);
        } else if(opts.kmer_result_ == FULL_SETSKETCH)
            make_save(fss);
    } else if(opts.sspace_ == SPACE_MULTISET) make_save(bmhs);
    else if(opts.sspace_ == SPACE_PSET) make(pmhs);
    else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        if(opts.parse_by_seq_) {
            omhs.reserve(nt);
            for(size_t i = 0; i < nt; omhs.emplace_back(ss, opts.k_), ++i);
        } else {
            throw std::invalid_argument("Space edit distance is only available in parse-by-seq mode, as it is only defined on strings rather than string collections.");
        }
    }
    while(ctrs.size() < nt) ctrs.emplace_back(opts.cssize());
#define __RESET(tid) do { \
        if(!opss.empty()) opss[tid].reset();\
        else if(!fss.empty()) fss[tid].reset();\
        else if(!bmhs.empty()) bmhs[tid].reset();\
        else if(!pmhs.empty()) pmhs[tid].reset();\
        /*else if(!omhs.empty()) omhs[tid].reset();*/\
        if(ctrs.size() > unsigned(tid)) ctrs[tid].reset();\
    } while(0)

    if(opts.parse_by_seq_) {
        std::vector<FastxSketchingResult> res(paths.size());
        ret.nperfile_.resize(paths.size());
        if(paths.size() == 1) {
            res.front() = fastx2sketch_byseq(opts, paths.front(), kseqs.kseqs_, true);
        } else {
            OMP_PFOR_DYN
            for(size_t i = 0; i < paths.size(); ++i) {
                std::fprintf(stderr, "sketching file %s at idx %zu\n", paths[i].data(), i);
                res[i] = fastx2sketch_byseq(opts, paths[i], kseqs.kseqs_);
                std::fprintf(stderr, "Sketched %zu/%zu (%s)\n", i + 1, paths.size(), paths[i].data());
            }
        }
        //std::fprintf(stderr, "Merging files\n");
        ret = FastxSketchingResult::merge(res.data(), res.size(), paths);
    } else {
        if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            throw std::runtime_error("edit distance is only available in parse by seq mode");
        }
        if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET) {
             opts.save_kmercounts_ = true; // Always save counts for PMinHash and BagMinHash
        }
        if(paths.size() == 1)
            std::fprintf(stderr, "Currently, only one thread is used per file in sketching. This may be slow if only one file is being processed.\n");
        ret.destination_files_.resize(paths.size());
        if(opts.save_kmers_) {
            ret.kmerfiles_.resize(paths.size());
        }
        if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
            ret.kmercountfiles_.resize(paths.size());
        }
        ret.cardinalities_.resize(paths.size());
#ifndef NDEBUG
        for(size_t i = 0; i < ret.names_.size(); ++i) {
            std::fprintf(stderr, "name %zu is %s\n", i, ret.names_[i].data());
        }
#endif
        std::fprintf(stderr, "kmer result type: %s\n", to_string(opts.kmer_result_).data());
        std::fprintf(stderr, "sketching space type: %s\n", to_string(opts.sspace_).data());
        std::string suffix = to_suffix(opts);
        auto makedest = [&](const std::string &path) -> std::string {
            std::string ret(path);
            ret = ret.substr(0, ret.find_first_of(' '));
            if(opts.trim_folder_paths_ || opts.outprefix_.size()) {
                ret = trim_folder(path);
                if(opts.outprefix_.size())
                    ret = opts.outprefix_ + '/' + ret;
            }
            if(opts.kmer_result_ <= FULL_SETSKETCH)
                ret = ret + std::string(".sketchsize") + std::to_string(opts.sketchsize_);
            ret = ret + std::string(".k") + std::to_string(opts.k_);
            if(opts.w_ > opts.k_) {
                ret = ret + std::string(".w") + std::to_string(opts.w_);
            }
            if(opts.count_threshold_ > 0) {
                ret = ret + ".ct_threshold";
                if(std::fmod(opts.count_threshold_, 1.)) ret = ret + std::to_string(opts.count_threshold_);
                else ret = ret + std::to_string(int(opts.count_threshold_));
            }
            if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_EDIT_DISTANCE) {
                ret += '.';
                ret += to_string(opts.ct());
                if(opts.ct() != EXACT_COUNTING)
                    ret += std::to_string(opts.cssize_);
            }
            if(opts.kmer_result_ <= FULL_SETSKETCH) {
                ret = ret + "." + to_string(opts.sspace_);
            } else {
                ret = ret + "." + to_string(opts.kmer_result_);
            }
            ret = ret + "." + bns::to_string(opts.rht_);
            ret = ret + suffix;
            return ret;
        };
        if(opts.build_sig_matrix_) {
            ret.signatures_.resize(ss * paths.size());
        }
        if((opts.build_mmer_matrix_ || opts.save_kmers_) && opts.kmer_result_ != FULL_MMER_SEQUENCE) {
            ret.kmers_.resize(ss * paths.size());
        }
        if(opts.build_count_matrix_) {
            ret.kmercounts_.resize(ss * paths.size());
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < paths.size(); ++i) {
            int tid = 0;
            OMP_ONLY(tid = omp_get_thread_num();)
            //const int tid = OMP_ELSE(omp_get_thread_num(), 0);
            //const auto starttime = std::chrono::high_resolution_clock::now();
            auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
            const size_t mss = ss * myind;
            auto &path = paths[myind];
            //std::fprintf(stderr, "parsing from path = %s\n", path.data());
            ret.destination_files_[myind] = makedest(path);
            auto &destination = ret.destination_files_[myind];
            const std::string destination_prefix = destination.substr(0, destination.find_last_of('.'));
            std::string destkmercounts = destination_prefix + ".kmercounts.f64";
            std::string destkmer = destination_prefix + ".kmer.u64";
            const bool dkif = check_compressed(destkmer);
            const bool dkcif = check_compressed(destkmercounts);
            const bool destisfile = check_compressed(destination);
            if(ret.kmercountfiles_.size() > myind) ret.kmercountfiles_[myind] = destkmercounts;
            if(opts.cache_sketches_ &&
               (destisfile || (opts.kmer_result_ == FULL_MMER_COUNTDICT && dkif)) &&
               (!opts.save_kmers_ || dkif) &&
               ((!opts.save_kmercounts_ && opts.kmer_result_ != FULL_MMER_COUNTDICT) || dkcif)
            )
            {
                if(opts.kmer_result_ < FULL_MMER_SET) {
                    if(ret.signatures_.size()) {
                        load_copy(destination, &ret.signatures_[mss]);
                        ret.cardinalities_[myind] = compute_cardest(&ret.signatures_[mss], ss);
                    }
                    if(ret.kmers_.size())
                        load_copy(destkmer, &ret.kmers_[mss]);
                    if(ret.kmercounts_.size())
                        load_copy(destkmercounts, &ret.kmercounts_[mss]);
                } else if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                    if(!dkcif)
                        throw std::runtime_error(std::string("Expected destkmercounts (") + destkmercounts + ") to be a file. Run again?");
                    //std::fprintf(stderr, "Cached at path %s, %s, %s\n", destination.data(), destkmercounts.data(), destkmer.data());
                    if(!iscomp(destkmercounts)) {
                        //std::fprintf(stderr, "Attempting to memory-map file at %s\n", destkmercounts.data());
                        mio::mmap_sink ms(destkmercounts);
                        if(ms.size() % sizeof(double)) throw std::runtime_error(std::string("Wrong size file ") + destkmercounts);
                        ret.cardinalities_[myind] = std::accumulate((const double *)ms.data(),(const double *)&ms[ms.size()], 0.);
                        continue;
                    }
                    const auto dkcsz = destkmercounts.size();
                    std::string cmd;
                    if(std::equal(&destkmercounts[dkcsz - 3], &destkmercounts[dkcsz], ".gz")) {
                        cmd = "gzip -dc ";
                    } else if(std::equal(&destkmercounts[dkcsz - 3], &destkmercounts[dkcsz], ".xz")) {
                        cmd = "xz -dc ";
                    } else if(destkmercounts.size() > 4 && std::equal(&destkmercounts[dkcsz - 4], &destkmercounts[dkcsz], ".bz2")) {
                        cmd = "bzip2 -dc ";
                    } else __builtin_unreachable();
                    cmd = cmd + destkmercounts;
                    std::fprintf(stderr, "Command '%s'\n", cmd.data());
                    std::FILE *ifp = ::popen(cmd.data(), "r");
                    if(!ifp) throw std::runtime_error("Failed to perform popen call for decompressing k-mer counts");
                    static constexpr size_t N = 4096;
                    size_t n = 4096;
                    std::unique_ptr<double[]> d(new double[N]);
                    long double s = 0.;
                    for(;n == N;) {
                        n = std::fread(d.get(), sizeof(double), N, ifp);
                        s = std::accumulate(d.get(), &d[n], s);
                    }
                    ::pclose(ifp);
                    ret.cardinalities_[myind] = s;
                } else if(opts.kmer_result_ == FULL_MMER_SET) {
                    ret.cardinalities_[myind] = bns::filesize(path.data()) / (opts.use128() ? 16: 8);
                }
                DBG_ONLY(std::fprintf(stderr, "Cache-sketches enabled. Using saved data at %s\n", destination.data());)
                continue;
            } else {
#ifndef NDEBUG
                std::fprintf(stderr, "We skipped caching because: %d is cache sketches\n", opts.cache_sketches_);
                std::fprintf(stderr, "destisfile: %d. is countdict %d. is kmerfile %d\n", destisfile, opts.kmer_result_ == FULL_MMER_COUNTDICT, dkif);
                std::fprintf(stderr, "kc save %d, kmer result %s, dkcif %d\n", opts.save_kmercounts_, to_string(opts.kmer_result_).data(), dkcif);
#endif
            }
            __RESET(tid);
            auto perf_for_substrs = [&](const auto &func) {
                for_each_substr([&](const std::string &subpath) {
                    //std::fprintf(stderr, "Doing for_each_substr for subpath = %s\n", subpath.data());
                    auto lfunc = [&](auto x) {if(!opts.fs_ || !opts.fs_->in_set(x)) func(x);};
#define FUNC_FE(f) f(lfunc, subpath.data(), kseqs.kseqs_ + tid)
                    if(!opts.parse_protein() && (opts.w_ > opts.k_ || opts.k_ <= 64)) {
                        if(opts.k_ < 32) {
                            auto encoder(opts.enc_);
                            FUNC_FE(encoder.for_each);
                        } else {
                            auto encoder(opts.enc_.to_u128());
                            FUNC_FE(encoder.for_each);
                        }
                    } else if(opts.use128()) {
                        DBG_ONLY(std::fprintf(stderr, "Parsing Protein with k = %u for 128-bit hashes\n", opts.k_);)
                        FUNC_FE(opts.rh128_.for_each_hash);
                    } else {
                        DBG_ONLY(std::fprintf(stderr, "Parsing Protein with k = %u for 64-bit hashes\n", opts.k_);)
                        FUNC_FE(opts.rh_.for_each_hash);
                    }
#undef FUNC_FE
                }, path);
            };
            const bool setsketch_with_counts = (opts.kmer_result_ == FULL_SETSKETCH) && (opts.save_kmercounts_ || opts.count_threshold_ > 0);
            if(
                (opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT)
                 || setsketch_with_counts
            )
            {
                auto &ctr = ctrs[tid];
                perf_for_substrs([&ctr](auto x) {ctr.add(x);});
                std::vector<u128_t> kmervec128;
                std::vector<uint64_t> kmervec64;
                std::vector<double> kmerveccounts;
                if(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                    if(opts.use128()) {
                        ctr.finalize(kmervec128, kmerveccounts, opts.count_threshold_);
                    } else {
                        ctr.finalize(kmervec64, kmerveccounts, opts.count_threshold_);
                    }
                    ret.cardinalities_[myind] = opts.kmer_result_ == FULL_MMER_SET ? (opts.use128() ? kmervec128.size(): kmervec64.size()): std::accumulate(kmerveccounts.begin(), kmerveccounts.end(), size_t(0));
                } else if(opts.sspace_ == SPACE_MULTISET) {
                    ctr.finalize(bmhs[tid], opts.count_threshold_);
                    ret.cardinalities_[myind] = bmhs[tid].total_weight();
                    std::copy(bmhs[tid].data(), bmhs[tid].data() + ss, &ret.signatures_[mss]);
                } else if(opts.sspace_ == SPACE_PSET) {
                    ctr.finalize(pmhs[tid], opts.count_threshold_);
                    std::copy(pmhs[tid].data(), pmhs[tid].data() + ss, &ret.signatures_[mss]);
                    ret.cardinalities_[myind] = pmhs[tid].total_weight();
                } else if(setsketch_with_counts) {
                    assert(fss.size());
                    ctr.finalize(fss[tid], opts.count_threshold_);
                    ret.cardinalities_[myind] = fss[tid].getcard();
                } else throw std::runtime_error("Unexpected space for counter-based m-mer encoding");
                    // Make bottom-k if we generated full k-mer sets or k-mer count dictionaries, and copy then over
                if(kmervec64.size() || kmervec128.size()) {
                    //std::fprintf(stderr, "If we gathered full k-mers, and we asked for signatures, let's store bottom-k hashes in the signature space\n");
                    if(ret.signatures_.size()) {
                        std::vector<BKRegT> keys(ss);
                        auto kvcp = kmerveccounts.data();
                        if(kmerveccounts.empty()) kvcp = nullptr;
                        if(kmervec128.size()) bottomk(kmervec128, keys, opts.count_threshold_, kvcp);
                        else bottomk(kmervec64, keys, opts.count_threshold_, kvcp);
                        std::copy(keys.begin(), keys.end(), (BKRegT *)&ret.signatures_[mss]);
                    }
                }
                std::FILE * ofp = std::fopen(destination.data(), "wb");
                if(!ofp) throw std::runtime_error(std::string("Failed to open std::FILE * at") + destination);
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
                } else if((opts.kmer_result_ ==  ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH)) {
                    buf = (const void *)fss[tid].data();
                    nb = ss * sizeof(RegT);
                    srcptr = fss[tid].data();
                } else nb = 0, srcptr = nullptr;
                if(srcptr && ret.signatures_.size())
                    std::copy(srcptr, srcptr + ss, &ret.signatures_[mss]);
                checked_fwrite(ofp, buf, nb);
                std::fclose(ofp);
                if((opts.save_kmers_ || opts.build_mmer_matrix_) && !(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_COUNTDICT)) {
                    assert(ret.kmerfiles_.size());
                    ret.kmerfiles_[myind] = destkmer;
                    const uint64_t *ptr = opts.sspace_ == SPACE_PSET ? pmhs[tid].ids().data():
                                      opts.sspace_ == SPACE_MULTISET ? bmhs[tid].ids().data():
                                      opts.kmer_result_ == ONE_PERM ? opss[tid].ids().data() :
                                      opts.kmer_result_ == FULL_SETSKETCH ? fss[tid].ids().data():
                                          static_cast<uint64_t *>(nullptr);
                    if(!ptr) throw 2;
                    if((ofp = std::fopen(destkmer.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer file");
                    //std::fprintf(stderr, "Writing to file %s\n", destkmer.data());

                    checked_fwrite(ofp, ptr, sizeof(uint64_t) * ss);
                    DBG_ONLY(std::fprintf(stderr, "About to copy to kmers of size %zu\n", ret.kmers_.size());)
                    if(ret.kmers_.size())
                        std::copy(ptr, ptr + ss, &ret.kmers_[mss]);
                    if(ofp) std::fclose(ofp);
                }
                if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                    //std::fprintf(stderr, "About to save kmer counts manually\n");
                    assert(ret.kmercountfiles_.size());
                    ret.kmercountfiles_.at(i) = destkmercounts;
                    if((ofp = std::fopen(destkmercounts.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer counts");
                    std::vector<double> tmp(ss);
#define DO_IF(x) if(x.size()) {std::copy(x[tid].idcounts().begin(), x[tid].idcounts().end(), tmp.data());}
                    if(opts.kmer_result_ == FULL_MMER_COUNTDICT || (opts.kmer_result_ == FULL_MMER_SET && opts.save_kmercounts_)) {
                        std::fprintf(stderr, "kvc size %zu. Writing to file %s\n", kmerveccounts.size(), destkmercounts.data());
                        tmp = std::move(kmerveccounts);
                    } else DO_IF(pmhs) else DO_IF(bmhs) else DO_IF(opss) else DO_IF(fss)
#undef DO_IF
                    const size_t nb = tmp.size() * sizeof(double);
                    checked_fwrite(ofp, tmp.data(), nb);
                    std::fclose(ofp);
                    if(ret.kmercounts_.size()) {
                        std::fprintf(stderr, "Copying range of size %zu from tmp to ret.kmercounts of size %zu\n", tmp.size(), ret.kmercounts_.size());
                        std::copy(tmp.begin(), tmp.begin() + ss, &ret.kmercounts_[mss]);
                    }
                }
            } else if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                ret.kmers_.clear();
                //std::fprintf(stderr, "Full mmer sequence\n");
                std::FILE * ofp;
                if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                void *dptr = nullptr;
                size_t m = 1 << 20;
                size_t l = 0;
                if(posix_memalign(&dptr, 16, (1 + opts.use128()) * m * sizeof(uint64_t))) throw std::bad_alloc();

                perf_for_substrs([&](auto x) {
                    using DT = decltype(x);
                    DT *ptr = (DT *)dptr;
                    if(opts.homopolymer_compress_minimizers_ && l > 0 && ptr[l - 1] == x) return;
                    if(l == m) {
                        size_t newm = m << 1;
                        void *newptr = nullptr;
                        if(posix_memalign((void **)&newptr, 16, newm * sizeof(DT))) throw std::bad_alloc();
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
                // These occur twice, because if the user asks for counts, or if the user asks for a minimum count level for inclusion.
                // Because of this, we have to generate the key-count map.
                // Those cases are handled above with the count-based methods.
                //std::fprintf(stderr, "kmer result is oneperm or setsketch. Dest: %s\n", destination.data());
                std::FILE * ofp;
                if((ofp = std::fopen(destination.data(), "wb")) == nullptr)
                    throw std::runtime_error(std::string("Failed to open file") + destination + "for writing minimizer sequence");
                if(opss.empty() && fss.empty()) throw std::runtime_error("Both opss and fss are empty\n");
                const size_t opsssz = opss.size();
                if(opsssz) {
                    //std::fprintf(stderr, "Encode for the opset sketch, %zu is the size for tid %d\n", opss.size(), tid);
                    assert(opss.size() > unsigned(tid));
                    assert(opss.at(tid).total_updates() == 0);
                    auto p = &opss[tid];
                    perf_for_substrs([p](auto hv) {p->update(hv);});
                    //std::fprintf(stderr, "Encode for the opset sketch. card now: %g, %zu updates\n", opss[tid].getcard(), opss[tid].total_updates());
                    assert(ret.cardinalities_.size() > i);
                    ret.cardinalities_[myind] = p->getcard();
                } else {
                    //std::fprintf(stderr, "Encode for the set sketch\n");
                    perf_for_substrs([p=&fss[tid]](auto hv) {p->update(hv);});
                    ret.cardinalities_[myind] = fss[tid].getcard();
                }
                const uint64_t *ids = nullptr;
                const uint32_t *counts = nullptr;
                const RegT *ptr = opsssz ? opss[tid].data(): fss[tid].data();
                assert(ptr);
                if(opts.build_mmer_matrix_)
                    ids = opsssz ? opss[tid].ids().data(): fss[tid].ids().data();
                if(opts.build_count_matrix_)
                    counts = opsssz ? opss[tid].idcounts().data(): fss[tid].idcounts().data();
                ::write(::fileno(ofp), ptr, sizeof(RegT) * ss);
                //checked_fwrite(ofp, ptr, sizeof(RegT) * ss);
                std::fclose(ofp);
                if(ptr && ret.signatures_.size()) std::copy(ptr, ptr + ss, &ret.signatures_[mss]);
                if(ids && ret.kmers_.size())
                    std::copy(ids, ids + ss, &ret.kmers_[mss]);
                if(counts && ret.kmercounts_.size())
                    std::copy(counts, counts + ss, &ret.kmercounts_[mss]);
            } else throw std::runtime_error("Unexpected: Not FULL_MMER_SEQUENCE, FULL_MMER_SET, ONE_PERM, FULL_SETSKETCH, SPACE_MULTISET, or SPACE_PSET");
        } // parallel paths loop
        ret.names_ = paths;
    }
    return ret;
}
SketchingResult SketchingResult::merge(SketchingResult *start, size_t n, const std::vector<std::string> &names={}) {
    std::fprintf(stderr, "About to merge from %p of size %zu, names has size %zu\n", (void *)start, n, names.size());
    SketchingResult ret;
    ret.options_ = start->options_;
    if(n == 0) return ret;
    else if(n == 1) {
        ret = std::move(*start);
        std::transform(ret.names_.begin(), ret.names_.end(), ret.names_.begin(), [&names](const auto &x) {return names.front() + ":" + x;});
        return ret;
    }
    //ret.nperfile_.resize(total_seq);
    for(size_t i = 0; i < n; ++i) {
        ret.nperfile_.insert(ret.nperfile_.end(), start[i].nperfile_.begin(), start[i].nperfile_.end());
    }
    size_t total_seqs = 0, total_sig_size = 0;
    std::vector<size_t> offsets(n + 1);
    std::vector<size_t> sig_offsets(n + 1);
    for(size_t i = 0; i < n; ++i) {
        const size_t nseqsi = start[i].names_.size();
        const size_t nregs = start[i].signatures_.size();
        assert(start[i].names_.size() == start[i].cardinalities_.size());
        total_seqs += nseqsi;
        total_sig_size += nregs;
        offsets[i + 1] = total_seqs;
        sig_offsets[i + 1] = total_sig_size;
    }
    ret.names_.resize(total_seqs);
    if(std::any_of(start, start + n, [](auto &x) {return x.sequences_.size();})) {
        ret.sequences_.resize(total_seqs);
    }
    const size_t sketchsz = start->signatures_.size() / start->names_.size();
    if(total_sig_size > 0) {
        ret.signatures_.resize(total_sig_size);
    }
    if(start->kmers_.size()) {
        ret.kmers_.resize(total_seqs * sketchsz);
    }
    ret.cardinalities_.resize(total_seqs);
    if(start->kmercounts_.size()) {
        ret.kmercounts_.resize(total_sig_size);
    }
    const bool seqsz = total_seqs,
               kmercountsz = !start->kmercounts_.empty();
    for(size_t i = 0; i < n; ++i) {
        auto &src = start[i];
        assert(src.names_.size() == offsets[i + 1] - offsets[i]);
        const auto ofs = offsets[i], sofs = sig_offsets[i];
        std::string fname;
        if(names.size() > i) fname = names[i].substr(0, names[i].find_first_of(' '));
        // Append filename to sequence names to ensure seq names
        std::transform(src.names_.begin(), src.names_.end(), &ret.names_[ofs], [&fname](const auto &x) {
            return x + ':' + fname;
        });
        //std::fprintf(stderr, "Copy cardinalities of size %zu to %zu/%zu\n", src.cardinalities_.size(), ofs, ret.cardinalities_.size());
        std::copy(src.cardinalities_.begin(), src.cardinalities_.end(), &ret.cardinalities_.at(ofs));
        //std::fprintf(stderr, "Copy sequences\n");
        if(seqsz) {
            //std::fprintf(stderr, "Copying %zd sequences to idx %zu/%zu\n", src.sequences_.size(), ofs, ret.sequences_.size());
            std::copy(src.sequences_.begin(), src.sequences_.end(), &ret.sequences_.at(ofs));
        }
        if(!start[i].signatures_.empty())
            std::copy(src.signatures_.begin(), src.signatures_.end(), &ret.signatures_.at(sofs));
        if(!start[i].kmers_.empty())
            std::copy(src.kmers_.begin(), src.kmers_.end(), &ret.kmers_[sofs]);
        if(kmercountsz)
            std::copy(src.kmercounts_.begin(), src.kmercounts_.end(), &ret.kmercounts_[sofs]);
    }
    return ret;
}



} // dashing2
