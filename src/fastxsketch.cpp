#include "fastxsketch.h"
#include "mio.hpp"
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
    for(size_t i = k; i > 0;) {
        --i;
        ret[i] = pq.top();
        pq.pop();
    }
}

template<typename T>
void load_copy(const std::string &path, T *ptr) {
    std::FILE *fp = std::fopen(path.data(), "rb");
    if(!fp) throw std::runtime_error(std::string("Failed to open ") + path);
    struct stat st;
    if(::fstat(::fileno(fp), &st)) {
        throw std::runtime_error(std::string("Failed to get fd from ") + path + std::to_string(::fileno(fp)));
    }
    std::fprintf(stderr, "Size of file at %s: %zu\n", path.data(), size_t(st.st_size));
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

FastxSketchingResult fastx2sketch(Dashing2Options &opts, const std::vector<std::string> &paths) {
    if(paths.empty()) throw std::invalid_argument("Can't sketch empty path set");
    const size_t nt = opts.nthreads();
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
        if(opts.kmer_result_ == ONE_PERM) make_save(opss);
        else if(opts.kmer_result_ == FULL_SETSKETCH) {
            make_save(fss);
        }
    } else if(opts.sspace_ == SPACE_MULTISET) {
        make_save(bmhs);
    } else if(opts.sspace_ == SPACE_PSET) {
        make(pmhs);
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        if(opts.parse_by_seq_) {
            omhs.reserve(nt);
            for(size_t i = 0; i < nt; omhs.emplace_back(ss, opts.k_), ++i);
        } else {
            throw std::invalid_argument("Space edit distance is only available in parse-by-seq mode, as it is only defined on strings rather than string collections.");
        }
    }
    while(ctrs.size() < nt) ctrs.emplace_back(opts.cssize());
    auto reset = [&](int tid) {
        if(!fss.empty()) fss[tid].reset();
        else if(!opss.empty()) opss[tid].reset();
        else if(!bmhs.empty()) bmhs[tid].reset();
        else if(!pmhs.empty()) pmhs[tid].reset();
        else if(!omhs.empty()) omhs[tid].reset();
        //else throw std::runtime_error("Unexpected: no sketches are available");
        if(ctrs.size() > unsigned(tid)) ctrs[tid].reset();
    };
    if(opts.parse_by_seq_) {
        std::vector<FastxSketchingResult> res(paths.size());
        ret.nperfile_.resize(paths.size());
        OMP_PFOR_DYN
        for(size_t i = 0; i < paths.size(); ++i) {
            res[i] = fastx2sketch_byseq(opts, paths[i], kseqs.kseqs_);
            DBG_ONLY(std::fprintf(stderr, "Sketched %zu/%zu (%s)\n", i + 1, paths.size(), paths[i].data());)
        }
        DBG_ONLY(std::fprintf(stderr, "Merging files\n");)
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
        ret.names_.resize(paths.size());
        ret.cardinalities_.resize(paths.size());
        std::copy(paths.begin(), paths.end(), ret.names_.begin());
        for(size_t i = 0; i < ret.names_.size(); ++i) {
            std::fprintf(stderr, "name %zu is %s\n", i, ret.names_[i].data());
        }
        std::fprintf(stderr, "kmer result type: %s\n", to_string(opts.kmer_result_).data());
        std::fprintf(stderr, "sketching space type: %s\n", to_string(opts.sspace_).data());
        std::string suffix = to_suffix(opts);
        auto makedest = [&](const std::string &path) -> std::string {
            std::string ret(path);
            ret = ret.substr(0, ret.find_first_of(' '));
            if(opts.trim_folder_paths_) {
                ret = trim_folder(path);
                if(opts.outprefix_.size()) {
                    ret = opts.outprefix_ + '/' + ret;
                }
            }
            ret = ret + std::string(".") + std::to_string(opts.k_);
            if(opts.w_ > opts.k_) {
                ret = ret + std::string(".") + std::to_string(opts.w_);
            }
            ret = ret + "." + bns::to_string(opts.rht_);
            ret = ret + suffix;
            return ret;
        };
        if(opts.build_sig_matrix_) {
            ret.signatures_.resize(ss * paths.size());
        }
        if(opts.build_mmer_matrix_ || opts.save_kmers_) {
            ret.kmers_.resize(ss * paths.size());
        }
        std::fprintf(stderr, "ret kmer size %zu\n", ret.kmers_.size());
        if(opts.build_count_matrix_) {
            ret.kmercounts_.resize(ss * paths.size());
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < paths.size(); ++i) {
            const int tid = OMP_ELSE(omp_get_thread_num(), 0);
            std::fprintf(stderr, "Sketching from tid %d at index %zu\n", tid, i);
            const auto starttime = std::chrono::high_resolution_clock::now();
            auto &path = paths[i];
            const std::string destination = makedest(path); // Unused right now, but soon won't be (?)
            const std::string destination_prefix = destination.substr(0, destination.find_last_of('.'));
            std::string destkmer = destination_prefix + ".kmer.u64";
            const std::string destkmercounts = destination_prefix + ".kmercounts.f64";
            const bool dkif = bns::isfile(destkmer);
            const bool dkcif = bns::isfile(destkmercounts);
            ret.destination_files_[i] = destination;
            if(ret.kmercountfiles_.size() > i) ret.kmercountfiles_[i] = destkmercounts;
            if(opts.cache_sketches_ &&
               bns::isfile(destination) &&
               (!opts.save_kmers_ || dkif) &&
               ((!opts.save_kmercounts_ && opts.kmer_result_ != FULL_MMER_COUNTDICT) || dkcif)
            )
            {
                if(opts.kmer_result_ < FULL_MMER_SET) {
                    if(ret.signatures_.size()) {
                        load_copy(destination, &ret.signatures_[ss * i]);
                        ret.cardinalities_[i] = compute_cardest(&ret.signatures_[ss * i], ss);
                    }
                    if(ret.kmers_.size())
                        load_copy(destkmer, &ret.kmers_[ss * i]);
                    if(ret.kmercounts_.size())
                        load_copy(destkmercounts, &ret.kmercounts_[ss * i]);
                } else if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                    assert(bns::isfile(destkmercounts));
                    mio::mmap_sink ms(destkmercounts);
                    const size_t nd = ms.size() / 8;
                    const double *p = (const double *)ms.data();
                    ret.cardinalities_[i] = std::accumulate(p, p + nd, 0.);
                    std::fprintf(stderr, "Setting card to %g from sum on disk from array of size %zu\n", ret.cardinalities_[i], nd);
                } else if(opts.kmer_result_ == FULL_MMER_SET) {
                    ret.cardinalities_[i] = bns::filesize(path.data());
                }
                std::fprintf(stderr, "Cache-sketches enabled. Using saved data at %s\n", destination.data());
                continue;
            }
            reset(tid);
            auto perf_for_substrs = [&](const auto &func) {
                for_each_substr([&](const std::string &subpath) {
                    auto lfunc = [&](auto x) {if(!opts.fs_ || !opts.fs_->in_set(x)) func(x);};
#define FUNC_FE(f) f(lfunc, subpath.data(), kseqs.kseqs_ + tid)
                    if(opts.w_ > opts.k_ || opts.k_ <= 64) {
                        if(opts.k_ < 32) {
                            auto encoder(opts.enc_);
                            FUNC_FE(encoder.for_each);
                        } else {
                            auto encoder(opts.enc_.to_u128());
                            FUNC_FE(encoder.for_each);
                        }
                    } else if(opts.use128()) {
                        auto hasher(opts.rh128_);
                        FUNC_FE(hasher.for_each_hash);
                    } else {
                        auto hasher(opts.rh_);
                        FUNC_FE(hasher.for_each_hash);
                    }
#undef FUNC_FE
                }, path);
            };
            const bool setsketch_with_counts = (opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) && (opts.save_kmercounts_ || opts.count_threshold_ > 0);
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
                    ret.cardinalities_[i] = opts.kmer_result_ == FULL_MMER_SET ? (opts.use128() ? kmervec128.size(): kmervec64.size()): std::accumulate(kmerveccounts.begin(), kmerveccounts.end(), size_t(0));
                } else if(opts.sspace_ == SPACE_MULTISET) {
                    ctr.finalize(bmhs[tid], opts.count_threshold_);
                    ret.cardinalities_[i] = bmhs[tid].total_weight();
                    std::copy(bmhs[tid].data(), bmhs[tid].data() + ss, &ret.signatures_[i * ss]);
                } else if(opts.sspace_ == SPACE_PSET) {
                    ctr.finalize(pmhs[tid], opts.count_threshold_);
                    std::copy(pmhs[tid].data(), pmhs[tid].data() + ss, &ret.signatures_[i * ss]);
                    ret.cardinalities_[i] = pmhs[tid].total_weight();
                } else if(setsketch_with_counts) {
                    if(opss.size())
                        ctr.finalize(opss[tid], opts.count_threshold_);
                    else
                        ctr.finalize(fss[tid], opts.count_threshold_);
                    ret.cardinalities_[i] = opss.size() ? opss[tid].getcard(): fss[tid].getcard();
                } else throw std::runtime_error("Unexpected space for counter-based m-mer encoding");
                    // Make bottom-k if we generated full k-mer sets or k-mer count dictionaries, and copy then over
                if(kmervec64.size() || kmervec128.size()) {
                    //std::fprintf(stderr, "If we gathered full k-mers, and we asked for signatures, let's store bottom-k hashes in the signature space\n");
                    if(ret.signatures_.size()) {
                        std::vector<BKRegT> keys(ss);
                        if(kmerveccounts.size()) {
                            auto kvp = kmerveccounts.data();
                            if(kmervec128.size()) bottomk(kmervec128, keys, opts.count_threshold_, kvp);
                            else bottomk(kmervec64, keys, opts.count_threshold_, kvp);
                        } else {
                            if(kmervec128.size()) bottomk(kmervec128, keys);
                            else bottomk(kmervec64, keys);
                        }
                        std::copy(keys.begin(), keys.end(), (BKRegT *)&ret.signatures_[i * ss]);
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
                    if(opss.size()) buf = (const void *)opss[tid].data();
                    else if(fss.size()) buf = (const void *)fss[tid].data();
                    else throw std::runtime_error("OPSS and FSS are both empty\n");
                    nb = ss * sizeof(RegT);
                    srcptr = opss.size() ? opss[tid].data(): fss[tid].data();
                } else nb = 0, srcptr = nullptr;
                if(srcptr && ret.signatures_.size())
                    std::copy(srcptr, srcptr + ss, &ret.signatures_[i * ss]);
                std::fprintf(stderr, "Copying out buffer of %zu to file %s\n", nb, destination.data());
                checked_fwrite(ofp, buf, nb);
                std::fclose(ofp);
                if((opts.save_kmers_ || opts.build_mmer_matrix_) && !(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_COUNTDICT)) {
                    assert(ret.kmerfiles_.size());
                    ret.kmerfiles_[i] = destkmer;
                    std::fprintf(stderr, "About to save kmers to %s\n", destkmer.data());
                    const uint64_t *ptr = opts.sspace_ == SPACE_PSET ? pmhs[tid].ids().data():
                                      opts.sspace_ == SPACE_MULTISET ? bmhs[tid].ids().data():
                                      opts.kmer_result_ == ONE_PERM ? opss[tid].ids().data() :
                                      opts.kmer_result_ == FULL_SETSKETCH ? fss[tid].ids().data():
                                          static_cast<uint64_t *>(nullptr);
                    if(!ptr) throw 2;
                    std::cerr << "Destkmer: " << destkmer << '\n';
                    if((ofp = std::fopen(destkmer.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer file");
                    std::fprintf(stderr, "Writing to file %s\n", destkmer.data());

                    checked_fwrite(ofp, ptr, sizeof(uint64_t) * ss);
                    std::fprintf(stderr, "About to copy to kmers of size %zu\n", ret.kmers_.size());
                    if(ret.kmers_.size())
                        std::copy(ptr, ptr + ss, &ret.kmers_[i * ss]);
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
                        tmp = kmerveccounts;
                        std::fprintf(stderr, "tmp size %zu, kvc size %zu. Writing to file %s\n", tmp.size(), kmerveccounts.size(), destkmercounts.data());
                    } else DO_IF(pmhs) else DO_IF(bmhs) else DO_IF(opss) else DO_IF(fss)
#undef DO_IF
                    const size_t nb = tmp.size() * sizeof(double);
                    checked_fwrite(ofp, tmp.data(), nb);
                    std::fclose(ofp);
                    if(ret.kmercounts_.size()) {
                        std::fprintf(stderr, "Copying range of size %zu from tmp to ret.kmercounts of size %zu\n", tmp.size(), ret.kmercounts_.size());
                        std::copy(tmp.begin(), tmp.begin() + ss, &ret.kmercounts_[i * ss]);
                    }
                }
            } else if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                std::fprintf(stderr, "Full mmer sequence\n");
                std::FILE * ofp;
                if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                // For faster build, we compile one loop.
                void *dptr = nullptr;
                size_t m = 1 << 20;
                size_t l = 0;
                if(posix_memalign((void **)&dptr, 16, (1 + opts.use128()) * m * sizeof(uint64_t))) throw std::bad_alloc();

                perf_for_substrs([&](auto x) {
                    using DT = decltype(x);
                    auto ptr = (DT *)dptr;
                    if(l == m) {
                        size_t newm = m << 1;
                        void *newptr = nullptr;
                        if(posix_memalign((void **)&newptr, 16, newm * sizeof(DT))) throw std::bad_alloc();
                        assert(newptr);
                        std::copy((DT *)dptr, (DT *)dptr + m, (DT *)newptr);
                        dptr = newptr;ptr = (DT *)dptr;
                        m = newm;
                    }
                    if(opts.homopolymer_compress_minimizers_) {
                        if(l > 0 && ptr[l - 1] == x) return;
                    }
                    ptr[l++] = x;
                });
                assert(dptr);
                checked_fwrite(ofp, dptr, l * (1 + opts.use128()) * sizeof(uint64_t));
                ret.cardinalities_[i] = l;
                std::free(dptr);
                std::fclose(ofp);
            } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
                // These occur twice, because if the user asks for counts, or if the user asks for a minimum count level for inclusion.
                // Because of this, we have to generate the key-count map.
                // Those cases are handled above with the count-based methods.
                //std::fprintf(stderr, "kmer result is oneperm or setsketch\n");
                std::FILE * ofp;
                if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                if(opss.empty() && fss.empty()) throw std::runtime_error("Both opss and fss are empty\n");
                const size_t opsssz = opss.size();
                perf_for_substrs([&](auto hv) {
                    if(opsssz) opss.operator[](tid).update(hv);
                    else fss.operator[](tid).update(hv);
                });
                const uint64_t *ids = nullptr;
                const uint32_t *counts = nullptr;
                const RegT *ptr = opsssz ? opss[tid].data(): fss[tid].data();
                ret.cardinalities_[i] = opsssz ? opss[tid].getcard(): fss[tid].getcard();
                assert(ptr);
                if(opts.build_mmer_matrix_)
                    ids = opsssz ? opss[tid].ids().data(): fss[tid].ids().data();
                if(opts.build_count_matrix_)
                    counts = opsssz ? opss[tid].idcounts().data(): fss[tid].idcounts().data();
                checked_fwrite(ofp, ptr, sizeof(RegT) * ss);
                std::fclose(ofp);
                if(ptr && ret.signatures_.size()) std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                if(ids && ret.kmers_.size())
                    std::copy(ids, ids + ss, &ret.kmers_[i * ss ]);
                if(counts && ret.kmercounts_.size())
                    std::copy(counts, counts + ss, &ret.kmercounts_[i * ss]);
            } else throw std::runtime_error("Unexpected: Not FULL_MMER_SEQUENCE, FULL_MMER_SET, ONE_PERM, FULL_SETSKETCH, SPACE_MULTISET, or SPACE_PSET");
            std::fprintf(stderr, "Sketching from tid %d at index %zu finished in %gms\n", tid, i, std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - starttime).count());;
        } // parallel paths loop
    }
    return ret;
}
SketchingResult SketchingResult::merge(SketchingResult *start, size_t n, const std::vector<std::string> &names={}) {
    SketchingResult ret;
    ret.options_ = start->options_;
    if(n == 0) return ret;
    else if(n == 1) return std::move(*start);
    ret.nperfile_.resize(n);
    size_t total_seqs = 0, total_sig_size = 0;
    std::vector<size_t> offsets(n + 1);
    std::vector<size_t> sig_offsets(n + 1);
    std::fprintf(stderr, "Computing offset arrays\n");
    for(size_t i = 0; i < n; ++i) {
        const size_t nseqsi = start[i].names_.size();
        const size_t nregs = start[i].signatures_.size();
        ret.nperfile_[i] = nseqsi;
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
               regsz = !start->signatures_.empty(),
               kmersz = !start->kmers_.empty(),
               kmercountsz = !start->kmercounts_.empty();
    std::fprintf(stderr, "Copying into merged thing with %zu total sequences\n", ret.names_.size());
    OMP_PFOR
    for(size_t i = 0; i < n; ++i) {
        auto &src = start[i];
        assert(src.names_.size() == offsets[i + 1] - offsets[i]);
        const auto ofs = offsets[i];
        std::string fname;
        if(names.size() > i) fname = names[i].substr(0, names[i].find_first_of(' '));
        // Append filename to sequence names to ensure seq names
        std::transform(src.names_.begin(), src.names_.end(), &ret.names_[ofs], [&fname](const auto &x) {
            return x + ':' + fname;
        });
        std::copy(src.cardinalities_.begin(), src.cardinalities_.end(), &ret.cardinalities_[ofs]);
        if(seqsz)
            std::copy(src.sequences_.begin(), src.sequences_.end(), &ret.sequences_[ofs]);
        if(regsz)
            std::copy(src.signatures_.begin(), src.signatures_.end(), &ret.signatures_[sig_offsets[i]]);
        if(kmersz)
            std::copy(src.kmers_.begin(), src.kmers_.end(), &ret.kmers_[sig_offsets[i]]);
        if(kmercountsz)
            std::copy(src.kmercounts_.begin(), src.kmercounts_.end(), &ret.kmercounts_[sig_offsets[i]]);
    }
    return ret;
}

struct OptSketcher {
    std::unique_ptr<BagMinHash> bmh;
    std::unique_ptr<ProbMinHash> pmh;
    std::unique_ptr<OPSetSketch> opss;
    std::unique_ptr<FullSetSketch> fss;
    std::unique_ptr<OrderMinHash> omh;
    std::unique_ptr<Counter> ctr;
    bns::RollingHasher<uint64_t> rh_;
    bns::RollingHasher<u128_t> rh128_;
    bns::Encoder<bns::score::Lex, uint64_t> enc_;
    bns::Encoder<bns::score::Lex, u128_t> enc128_;
    int k_, w_;
    bool use128_;
    OptSketcher(bns::Spacer sp, bool use128, bool canon): enc_(sp, canon), enc128_(std::move(enc_.to_u128())), k_(sp.k()), w_(sp.w()), use128_(use128) {
    }
    template<typename Func>
    void for_each(const Func &func, const char *s, size_t n) {
        if(k_ > 64) {
            if(use128_)
                rh128_.for_each_hash(func, s, n);
            else
                rh_.for_each_hash(func, s, n);
        } else if(k_ <= 32) {
            enc_.for_each(func, s, n);
        } else {
            enc128_.for_each(func, s, n);
        }
    }
    void reset() {
        if(bmh) bmh->reset();
        if(pmh) pmh->reset();
        if(opss) opss->reset();
        if(fss) fss->reset();
        if(omh) omh->reset();
        if(ctr) ctr->reset();
        rh_.reset();
        rh128_.reset();
    }
};

void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, std::vector<std::pair<char *, size_t>> &seqpairs, OptSketcher &sketchers) {
    const size_t oldsz = ret.names_.size();
    std::fprintf(stderr, "oldsize %zu, increment %zu\n", oldsz, newsz);
    newsz = oldsz + newsz;
    if(opts.build_sig_matrix_) {
        ret.signatures_.resize(newsz * opts.sketchsize_);
    }
    ret.cardinalities_.resize(newsz);
    std::fprintf(stderr, "mmer matrix size %zu. buuild %d, save kmers %d\n", ret.kmers_.size(), opts.build_mmer_matrix_, opts.save_kmers_);
    if(opts.build_mmer_matrix_ || opts.save_kmers_) {
        ret.kmers_.resize(opts.sketchsize_ * newsz);
    }
    if(opts.build_count_matrix_ || opts.save_kmercounts_) {
        ret.kmercounts_.resize(opts.sketchsize_ * newsz);
    }
    std::fprintf(stderr, "About to go through list, %zu sigsize, %zu mmer size %zu mmer count size\n", ret.signatures_.size(), ret.kmers_.size(), ret.kmercounts_.size());
    for(size_t i = oldsz; i < newsz; ++i) {
        std::fprintf(stderr, "%zu/%zu\n", i, newsz);
        auto seqp = seqpairs[i - oldsz];
        //std::fprintf(stderr, "Seq %zu/%s\n", seqp.second, seqp.first);
        sketchers.reset();
        if(sketchers.omh) {
            std::vector<uint64_t> res = sketchers.omh->hash(seqp.first, seqp.second);
            auto destp = &ret.signatures_[(i - oldsz) * opts.sketchsize_];
            if(sizeof(RegT) == sizeof(uint64_t)) { // double
                // Copy as-is
                std::memcpy(destp, res.data(), res.size() * sizeof(uint64_t));
            } else if(sizeof(RegT) > sizeof(uint64_t)) { // long double
                std::transform(res.begin(), res.end(), (u128_t *)destp, [](uint64_t x) {
                    long double ret;
                    std::memcpy(&ret, &x, sizeof(x));
                    ((uint64_t *)&ret)[1] = wy::wyhash64_stateless(&x);
                    return ret;
                });
            } else {
                // truncate to u32 from 64-bit hashes
                // these hashes are stronger than bottom-k because of positional information
                std::copy(res.begin(), res.end(), (uint32_t *)destp);
            }
        } else {
            sketchers.for_each([&](auto x) {
                if(opts.fs_ && opts.fs_->in_set(x)) return;
                if(sketchers.ctr) {
                    sketchers.ctr->add(x);
                } else if(sketchers.opss) {
                    sketchers.opss->update(x);
                } else if(sketchers.fss) {
                    sketchers.fss->update(x);
                } else throw std::runtime_error("Expected ctr or opss or fss and none are non-null");
            }, seqp.first, seqp.second);
            RegT *ptr = nullptr;
            const uint64_t *kmer_ptr = nullptr;
            std::vector<double> kmercounts;
            if(opts.sspace_ == SPACE_SET) {
                if(sketchers.ctr) {
                    double card;
                    if(sketchers.opss) {
                        sketchers.ctr->finalize(*sketchers.opss, opts.count_threshold_);
                        ptr = sketchers.opss->data();
                        card = sketchers.opss->getcard();
                    } else {
                        sketchers.ctr->finalize(*sketchers.fss, opts.count_threshold_);
                        ptr = sketchers.fss->data();
                        card = sketchers.fss->getcard();
                    }
                    ret.cardinalities_[i] = card;
                } else if(sketchers.opss) {
                    ptr = sketchers.opss->data();
                    ret.cardinalities_[i] = sketchers.opss->getcard();
                    if(sketchers.opss->ids().size())
                        kmer_ptr = sketchers.opss->ids().data();
                    if(sketchers.opss->idcounts().size()) {
                        kmercounts.resize(opts.sketchsize_);
                        std::copy(sketchers.opss->idcounts().begin(), sketchers.opss->idcounts().end(), kmercounts.begin());
                    }
                } else {
                    assert(sketchers.fss);
                    ptr = sketchers.fss->data();
                    ret.cardinalities_[i] = sketchers.fss->getcard();
                    if(sketchers.fss->ids().size())
                        kmer_ptr = sketchers.fss->ids().data();
                    if(sketchers.fss->idcounts().size()) {
                        kmercounts.resize(opts.sketchsize_);
                        std::copy(sketchers.fss->idcounts().begin(), sketchers.fss->idcounts().end(), kmercounts.begin());
                    }
                }
            } else if(opts.sspace_ == SPACE_MULTISET) {
                sketchers.ctr->finalize(*sketchers.bmh, opts.count_threshold_);
                ptr = sketchers.bmh->data();
                ret.cardinalities_[i] = sketchers.bmh->total_weight();
                if(sketchers.bmh->ids().size()) kmer_ptr = sketchers.bmh->ids().data();
                if(sketchers.bmh->idcounts().size()) {
                    kmercounts.resize(opts.sketchsize_);
                    std::copy(sketchers.bmh->idcounts().begin(), sketchers.bmh->idcounts().end(), kmercounts.begin());
                }
            } else if(opts.sspace_ == SPACE_PSET) {
                sketchers.ctr->finalize(*sketchers.pmh, opts.count_threshold_);
                ptr = sketchers.pmh->data();
                if(sketchers.pmh->ids().size()) kmer_ptr = sketchers.pmh->ids().data();
                ret.cardinalities_[i] = sketchers.pmh->total_weight();
                if(sketchers.pmh->idcounts().size()) {
                    kmercounts.resize(opts.sketchsize_);
                    std::copy(sketchers.pmh->idcounts().begin(), sketchers.pmh->idcounts().end(), kmercounts.begin());
                }
            } else throw std::runtime_error("Not yet implemented?");
            //std::fprintf(stderr, "Copying from %p to container of size %zu at %zu\n", (void *)ptr, ret.signatures_.size(), size_t(i * opts.sketchsize_));
            std::copy(ptr, ptr + opts.sketchsize_, &ret.signatures_[i * opts.sketchsize_]);
            if(kmer_ptr && ret.kmers_.size()) {
                std::copy(kmer_ptr, kmer_ptr + opts.sketchsize_, &ret.kmers_[i * opts.sketchsize_]);
            }
            if(kmercounts.size()) {
                std::copy(kmercounts.begin(), kmercounts.end(), &ret.kmercounts_[i * opts.sketchsize_]);
            }
        }
    }
    if(ret.sequences_.empty()) {
        for(auto &pair: seqpairs) std::free(pair.first);
    }
}

FastxSketchingResult fastx2sketch_byseq(Dashing2Options &opts, const std::string &path, kseq_t *kseqs, const size_t seqs_per_batch) {
    FastxSketchingResult ret;
    gzFile ifp;
    std::vector<std::pair<char *, size_t>> seqpairs(seqs_per_batch);
    kseq_t *myseq = nullptr;
    if(kseqs) myseq = &kseqs[OMP_ELSE(omp_get_thread_num(), 0)];
    else myseq = (kseq_t *)std::calloc(sizeof(kseq_t), 1);
    size_t batch_index = 0;
    OptSketcher sketching_data(opts.sp_, opts.use128(), opts.canonicalize());
    sketching_data.rh_ = opts.rh_;
    sketching_data.rh128_ = opts.rh128_;
    const bool save_ids = opts.save_kmers_ || opts.build_mmer_matrix_;
    const bool save_idcounts = opts.save_kmercounts_ || opts.build_count_matrix_;
    //std::fprintf(stderr, "save ids: %d, save counts %d\n", save_ids, save_idcounts);
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm_) sketching_data.opss.reset(new OPSetSketch(opts.sketchsize_, save_ids, save_idcounts));
        else sketching_data.fss.reset(new FullSetSketch(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_MULTISET) {
        sketching_data.bmh.reset(new BagMinHash(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_PSET) {
        sketching_data.pmh.reset(new ProbMinHash(opts.sketchsize_));
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        sketching_data.omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
    }
    if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || opts.count_threshold_ > 0.) {
        sketching_data.ctr.reset(new Counter(opts.cssize()));
    }

    for_each_substr([&](const auto &x) {
        if((ifp = gzopen(x.data(), "rb")) == nullptr) throw std::runtime_error(std::string("Failed to read from ") + x);
        bns::kseq_assign(myseq, ifp);
        for(int c;(c = kseq_read(myseq)) >= 0;) {
            ret.sequences_.emplace_back(myseq->seq.s, myseq->seq.l);
            char *seqstr = ret.sequences_.emplace_back(myseq->seq.s, myseq->seq.l).data();
            seqpairs[batch_index] = {seqstr, myseq->seq.l};
            ret.names_.emplace_back(myseq->name.s + (myseq->name.s[0] == '>'));
            if(++batch_index == seqs_per_batch) {
                resize_fill(opts, ret, batch_index, seqpairs, sketching_data);
                batch_index = 0;
            }
        }
        gzclose(ifp);
    }, path);
    if(batch_index) resize_fill(opts, ret, batch_index, seqpairs, sketching_data);
    if(!kseqs) kseq_destroy(myseq);
    return ret;
}
} // dashing2
