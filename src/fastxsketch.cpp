#include "d2.h"
#include <optional>
namespace dashing2 {

struct FastxSketchingResult {
    std::vector<std::string> names_;
    std::vector<std::string> destination_files_; // Contains sketches/kmer-sets,kmer-sequences etc.
    std::vector<std::string> kmerfiles_;         // Contains file-paths for k-mers, if saved.
    std::vector<std::string> kmercountfiles_;    // Contains k-mer counts, if saved
    std::vector<uint32_t> nperfile_; // This is either empty (in which case each filename/row has its own sketch)
                                     // Or, this contains a list indicating the number of sketches created for each file/line
    std::vector<RegT> signatures_; // Signatures, packed into a single array
    std::vector<uint64_t> kmers_;  // This contains the k-mers corresponding to signatures, if asked for
    std::vector<uint32_t> kmercounts_; // Contains counts for k-mers, if desired
};

struct SequenceSketch {
    // If counts were generated, rather than final sketches
    std::vector<Counter> ctrs;
    std::vector<RegT> sigs;
    std::vector<std::string> names;
};

template<typename F>
void for_each_substr(const F &func, const std::string &s, const int sep=' ') {
    const char *p;
    if((p = std::strchr(s.data(), sep)) == nullptr) {
        func(s);
        return;
    }
    const char *p2 = s.data();
    std::string tmp(p2, p);
    for(;;) {
        func(tmp);
        std::swap(p2, ++p);
        if((p = std::strchr(p2, sep)) == nullptr) {
            func(p2);
            break;
        }
        tmp = std::string(p2, p);
        if(std::all_of(tmp.begin(), tmp.end(), [](auto x) {return std::isspace(x);})) break;
    }
}

void checked_fwrite(std::FILE *fp, const void *ptr, size_t nb) {
    auto rc = std::fwrite(ptr, 1, nb, fp);
    if(rc != nb) throw std::runtime_error("Failed to buffered-write\n");
}

FastxSketchingResult fastx2sketch(ParseOptions &opts, std::vector<std::string> &paths) {
    if(paths.empty()) throw std::invalid_argument("Can't sketch empty path set");
    const size_t nt = opts.nthreads();
    const size_t ss = opts.sketchsize();
    FastxSketchingResult ret;
    std::vector<BagMinHash> bmhs;
    std::vector<ProbMinHash> pmhs;
    std::vector<OPSetSketch> opss;
    std::vector<FullSetSketch> fss;
    std::vector<OrderMinHash> omhs;
    std::vector<Counter> ctrs;
    kseq_t *kseqs = static_cast<kseq_t *>(std::calloc(nt, sizeof(kseq_t)));
    if(!kseqs) throw std::bad_alloc();
    auto make = [&](auto &x) {
        x.reserve(nt);
        for(size_t i = 0; i < nt; ++i)
            x.emplace_back(ss);
    };
    if(opts.sspace_ == SPACE_SET) {
        auto make_save = [&](auto &x) {
            x.reserve(nt);
            for(size_t i = 0; i < nt; ++i)
                x.emplace_back(ss, opts.save_kmers_, opts.save_kmercounts_);
        };
        if(opts.kmer_result_ == ONE_PERM) make_save(opss);
        else if(opts.kmer_result_ == FULL_SETSKETCH) {
            make_save(fss);
        }
    } else if(opts.sspace_ == SPACE_MULTISET) {
        make(bmhs);
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
        else throw std::runtime_error("Unexpected: no sketches are available");
    };
    if(opts.parse_by_seq_) {
        ret.nperfile_.resize(paths.size());
        std::fprintf(stderr, "Parse by seq not available; Dashing2 is still in development\n");
        std::exit(1);
    } else {
        if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            throw std::runtime_error("edit distance is only available in parse by seq mode");
        }
        if(paths.size() == 1) {
            // Parallelized construction
            throw std::runtime_error("Not yet completed: parallel sketching of a single file");
        } else {
            std::fprintf(stderr, "Multiple files, parsing each file/list of files with single thread. Total nthreads: %zu\n", nt);
            ret.destination_files_.resize(paths.size());
            if(opts.save_kmers_) {
                ret.kmerfiles_.resize(paths.size());
            }
            if(opts.save_kmercounts_) {
                ret.kmercountfiles_.resize(paths.size());
            }
            // TODO:
            // Change ordering to process larger datasets first for better load-balancing
            ret.names_.resize(paths.size());
            std::copy(paths.begin(), paths.end(), ret.names_.begin());
            std::string suffix = (opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) ?
                (opts.sspace_ == SPACE_SET ? ".ss": opts.sspace_ == SPACE_MULTISET ? ".bmh": opts.sspace_ == SPACE_PSET ? ".pmh": ".unknown")
                : (opts.kmer_result_ == FULL_MMER_SET ? ".kmerset" : opts.kmer_result_ == FULL_MMER_SEQUENCE ? ".mmerseq": ".unknown_kmer");
            if(opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_SET) {
                if(opts.use128())
                    suffix = suffix + "128";
                else
                    suffix = suffix + "64";
            }
            auto makedest = [&](const std::string &path) -> std::string {
                std::string ret(path);
                if(opts.trim_folder_paths_) {
                    ret = path.substr(0, path.find_last_of('/'));
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
            if(opts.build_mmer_matrix_) {
                ret.kmers_.resize(ss * paths.size());
            }
            if(opts.save_kmercounts_) {
                ret.kmercounts_.resize(ss * paths.size());
            }
            OMP_PFOR_DYN
            for(size_t i = 0; i < paths.size(); ++i) {
                const int tid = OMP_ELSE(omp_get_thread_num(), 0);
                kseq_t *const ksp = kseqs + tid;
                auto &path = paths[i];
                const std::string destination = makedest(path); // Unused right now, but soon won't be (?)
                const std::string destination_prefix = destination.substr(0, destination.find_last_of('.'));
                std::string destkmer = destination_prefix + ".kmer.u64";
                std::string destkmercounts = destination_prefix + ".kmercounts.f32";
                ret.destination_files_[i] = destination;
                reset(tid);
                if(opts.cache_sketches_ &&
                   bns::isfile(destination) &&
                   (!opts.save_kmers_ || bns::isfile(destkmer)) &&
                   (!opts.save_kmercounts_ || bns::isfile(destkmercounts))
                )
                {
                    std::fprintf(stderr, "Cache-sketches enabled. Using saved data at %s\n", destination.data());
                    continue;
                }
                reset(tid);
                auto perf_for_substrs = [&](const auto &func) {
                    for_each_substr([&](const std::string &subpath) {
                        if(opts.use128()) {
                            auto hasher(opts.rh128_);
                            hasher.for_each_hash(func, subpath.data(), ksp);
                        } else {
                            auto hasher(opts.rh_);
                            hasher.for_each_hash(func, subpath.data(), ksp);
                        }
                    }, path);
                };
                if(
                    (opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT)
                     ||
                    ((opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) && (opts.save_kmercounts_ || opts.count_threshold_ > 0))
                )
                {
                    auto &ctr = ctrs[tid];
                    perf_for_substrs([&ctr](auto x) {ctr.add(x);});
                    std::vector<u128_t> kmervec128;
                    std::vector<uint64_t> kmervec64;
                    std::vector<double> kmerveccounts;
                    if(opts.kmer_result_ == FULL_MMER_SET) {
                        if(opts.use128())
                            ctr.finalize(kmervec128, opts.count_threshold_);
                        else
                            ctr.finalize(kmervec64, opts.count_threshold_);
                    } else if(opts.sspace_ == SPACE_MULTISET) {
                        ctr.finalize(bmhs[tid], opts.count_threshold_);
                        auto ptr = bmhs[tid].data();
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    } else if(opts.sspace_ == SPACE_PSET) {
                        ctr.finalize(pmhs[tid], opts.count_threshold_);
                        auto ptr = pmhs[tid].data();
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    } else if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                        if(opts.use128())
                            ctr.finalize(kmervec128, kmerveccounts, opts.count_threshold_);
                        else
                            ctr.finalize(kmervec64, kmerveccounts, opts.count_threshold_);
                    } else throw std::runtime_error("Unexpected space for counter-based m-mer encoding");
                    std::FILE * ofp = std::fopen(destination.data(), "wb");
                    if(!ofp) throw std::runtime_error(std::string("Failed to open std::FILE * at") + destination);
                    const void *buf = nullptr;
                    long long int nb;
                    const RegT *srcptr;
                    if(kmervec128.size()) {
                        buf = (const void *)kmervec128.data();
                        nb = kmervec128.size() * sizeof(u128_t);
                    } else if(kmervec64.size()) {
                        buf = (const void *)kmervec64.data();
                        nb = kmervec64.size() * sizeof(uint64_t);
                    } else if(opts.sspace_ == SPACE_MULTISET) {
                        buf = (const void *)bmhs[tid].data();
                        nb = opts.sketchsize_ * sizeof(RegT);
                        srcptr = bmhs[tid].data();
                    } else if(opts.sspace_ == SPACE_PSET) {
                        buf = (const void *)pmhs[tid].data();
                        nb = opts.sketchsize_ * sizeof(RegT);
                        srcptr = pmhs[tid].data();
                    } else if((opts.kmer_result_ ==  ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH)) {
                        if(opss.size()) buf = (const void *)opss[tid].data();
                        else buf = (const void *)fss[tid].data();
                        nb = opts.sketchsize_ * sizeof(RegT);
                        srcptr = opss.size() ? opss[tid].data(): fss[tid].data();
                    } else nb = 0, srcptr = nullptr;
                    std::copy(srcptr, srcptr + ss, &ret.signatures_[i * ss]);
                    checked_fwrite(ofp, buf, nb);
                    std::fclose(ofp);
                    if(opts.save_kmers_) {
                        if(opts.sspace_ == SPACE_MULTISET)
                            throw std::runtime_error("Not yet available: tracking k-mer id's for BagMinHash");
                        assert(ret.kmerfiles_.size());
                        ret.kmerfiles_[i] = destkmer;
                        if((ofp = std::fopen(destkmer.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer file");

                        if(opts.sspace_ == SPACE_PSET) {
                            checked_fwrite(ofp, pmhs[tid].res_.data(), sizeof(pmhs[tid].res_[0]) * pmhs[tid].size());
                        } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
                            if(opss.size()) checked_fwrite(ofp, opss[tid].ids().data(), sizeof(opss[tid].ids()[0]) * opss[tid].ids().size());
                            else checked_fwrite(ofp, fss[tid].ids().data(), sizeof(fss[tid].ids()[0]) * fss[tid].ids().size());
                        }
                        if(ofp) std::fclose(ofp);
                    }
                    if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                        assert(ret.kmercountfiles_.size());
                        ret.kmercountfiles_[i] = destkmercounts;
                        if((ofp = std::fopen(destkmercounts.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer counts");
                        std::vector<float> tmp(opts.sketchsize_);
                        if(opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                            tmp.resize(kmerveccounts.size());
                            std::copy(kmerveccounts.data(), kmerveccounts.data() + kmerveccounts.size(), tmp.data());
                        } else if(pmhs.size()) {
                            auto p = pmhs[tid].idcounts().data();
                            std::copy(p, p + pmhs[tid].idcounts().size(), tmp.data());
                        } else if(bmhs.size()) {
                            auto p = pmhs[tid].idcounts().data();
                            std::copy(p, p + pmhs[tid].idcounts().size(), tmp.data());
                        } else if(opss.size()) {
                            auto p = opss[tid].idcounts().data();
                            std::copy(p, p + opss[tid].idcounts().size(), tmp.data());
                        } else if(fss.size()) {
                            auto p = fss[tid].idcounts().data();
                            std::copy(p, p + fss[tid].idcounts().size(), tmp.data());
                        }
                        checked_fwrite(ofp, tmp.data(), tmp.size() * sizeof(float));
                        std::fclose(ofp);
                    }
                } else if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
                    std::FILE * ofp;
                    if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                    // For faster compilation, we compile one loop
                    void *dptr;
                    size_t m = 1 << 20;
                    size_t l = 0;
                    if(posix_memalign((void **)&dptr, 16, (1 + opts.use128()) * m * sizeof(uint64_t))) throw std::bad_alloc();

                    perf_for_substrs([&](auto x) {
                        using DT = decltype(x);
                        auto ptr = (DT *)dptr;
                        if(l == m) {
                            void *newptr;
                            size_t newm = m << 1;
                            if(posix_memalign((void **)&newptr, 16, newm * sizeof(DT))) throw std::bad_alloc();
                            std::memcpy(newptr, dptr, (1 + opts.use128()) * m * sizeof(DT));
                            dptr = newptr;
                            m = newm;
                        }
                        if(opts.homopolymer_compress_minimizers_ && (l > 0 && ptr[l - 1] == x)) return;
                        ptr[l++] = x;
                    });
                    checked_fwrite(ofp, dptr, l * (1 + opts.use128()) * sizeof(uint64_t));
                    std::free(dptr);
                    std::fclose(ofp);
                } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
                    // These occur twice, because if the user asks for counts, or if the user asks for a minimum count level for inclusion.
                    // Because of this, we have to generate the key-count map.
                    // Those cases are handled above with the count-based methods.
                    const RegT *ptr = nullptr;
                    const uint64_t *ids = nullptr;
                    const uint32_t *counts = nullptr;
                    std::FILE * ofp;
                    if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                    perf_for_substrs([&](auto hv) {
                        if(opts.one_perm_) opss[tid].update(hv);
                        else fss[tid].update(hv);
                    });
                    ptr = opts.one_perm_ ? opss[tid].data(): fss[tid].data();
                    if(opts.save_kmers_ && opts.build_mmer_matrix_)
                        ids = opts.one_perm_ ? opss[tid].ids().data(): fss[tid].ids().data();
                    if(opts.save_kmercounts_) {
                        counts = opts.one_perm_ ? opss[tid].idcounts().data(): fss[tid].idcounts().data();
                    }
                    checked_fwrite(ofp, ptr, sizeof(RegT) * opts.sketchsize_);
                    std::fclose(ofp);
                    if(ret.signatures_.size()) std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    if(ids)
                        std::copy(ids, ids + ss, &ret.kmers_[i * ss ]);
                    if(counts)
                        std::copy(counts, counts + ss, &ret.kmercounts_[i * ss]);
                }
            }
        }
    }
    for(size_t i = 0; i < nt; ++i) bns::kseq_destroy_stack(kseqs[i]);
    std::free(kseqs);
    return ret;
}

} // dashing2
