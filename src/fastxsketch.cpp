#include "fastxsketch.h"
//#include <optional>
namespace dashing2 {


void checked_fwrite(std::FILE *fp, const void *ptr, size_t nb) {
    auto rc = std::fwrite(ptr, 1, nb, fp);
    if(rc != nb) throw std::runtime_error("Failed to buffered-write\n");
}

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
        for(const auto v: kmercounts_) {
            s += v; ss += v * double(v);
        }
        msg += "mean: ";
        msg += to_string(double(s / kcsz));
        std::cerr << msg << '\n';
        msg = msg + ", std " + to_string(double(std::sqrt(ss / kcsz - std::pow(s / kcsz, 2.))));
        std::cerr << msg << '\n';
    }
    return msg;
}

FastxSketchingResult fastx2sketch(ParseOptions &opts, std::vector<std::string> &paths) {
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
    kseq_t *kseqs = static_cast<kseq_t *>(std::calloc(nt, sizeof(kseq_t)));
    if(!kseqs) throw std::bad_alloc();
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
        ret.nperfile_.resize(paths.size());
        std::fprintf(stderr, "Parse by seq not available; Dashing2 is still in development\n");
        std::exit(1);
    } else {
        if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            throw std::runtime_error("edit distance is only available in parse by seq mode");
        }
        if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET) {
             opts.save_kmercounts_ = true; // Always save counts for PMinHash and BagMinHash
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
            if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                ret.kmercountfiles_.resize(paths.size());
            }
            // TODO:
            // Change ordering to process larger datasets first for better load-balancing
            ret.names_.resize(paths.size());
            std::copy(paths.begin(), paths.end(), ret.names_.begin());
            std::fprintf(stderr, "kmer result type: %s\n", to_string(opts.kmer_result_).data());
            std::string suffix = (opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) ?
                (opts.sspace_ == SPACE_SET ? (opts.one_perm_ ? ".opss": ".fss"):  opts.sspace_ == SPACE_MULTISET ? ".bmh": opts.sspace_ == SPACE_PSET ? ".pmh": ".unknown")
                : (opts.kmer_result_ == FULL_MMER_SET ? ".kmerset" : (opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_COUNTDICT) ? ".mmerseq": ".unknown_kmer");
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
            if(opts.build_count_matrix_) {
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
                std::fprintf(stderr, "Resetting\n");
                reset(tid);
                std::fprintf(stderr, "Resetting finished\n");
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
                const bool setsketch_with_counts = (opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) && (opts.save_kmercounts_ || opts.count_threshold_ > 0);
                if(
                    (opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT)
                     || setsketch_with_counts
                )
                {
                    std::fprintf(stderr, "Performing k-mer count\n");
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
                    } else if(opts.sspace_ == SPACE_MULTISET) {
                        ctr.finalize(bmhs[tid], opts.count_threshold_);
                        auto ptr = bmhs[tid].data();
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    } else if(opts.sspace_ == SPACE_PSET) {
                        ctr.finalize(pmhs[tid], opts.count_threshold_);
                        auto ptr = pmhs[tid].data();
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    } else if(setsketch_with_counts) {
                        std::fprintf(stderr, "Building setsketch with count threshold %g\n", opts.count_threshold_);
                        if(opts.kmer_result_ == ONE_PERM) {
                            ctr.finalize(opss[tid], opts.count_threshold_);
                        } else {
                            ctr.finalize(fss[tid], opts.count_threshold_);
                        }
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
                    long long int nb;
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
                    if(srcptr && ret.signatures_.size()) {
                        std::copy(srcptr, srcptr + ss, &ret.signatures_[i * ss]);
                    }
                    checked_fwrite(ofp, buf, nb);
                    std::fclose(ofp);
                    std::fprintf(stderr, "Save kmers\n");
                    if((opts.save_kmers_ || opts.build_mmer_matrix_) && !(opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_COUNTDICT)) {
                        assert(ret.kmerfiles_.size());
                        ret.kmerfiles_[i] = destkmer;
                        std::fprintf(stderr, "About to save kmers to %s\n", destkmer.data());
                        if((ofp = std::fopen(destkmer.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer file");

                        const uint64_t *ptr = nullptr;
                        if(opts.sspace_ == SPACE_PSET) {
                            static_assert(sizeof(pmhs[tid].res_[0]) == sizeof(uint64_t), "Must be 64-bit");
                            std::fprintf(stderr, "ProbMinHash\n");
                            ptr = pmhs[tid].ids().data();
                        } else if(opts.sspace_ == SPACE_MULTISET) {
                            std::fprintf(stderr, "Multiset sketch\n");
                            static_assert(sizeof(bmhs[tid].track_ids_[0]) == sizeof(uint64_t), "Must be 64-bit");
                            assert(bmhs[tid].ids().size());
                            ptr = bmhs[tid].ids().data();
                        } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
                            std::fprintf(stderr, "Sketch sketching. OPsz, Fsz %zu, %zu\n", opss.size(), fss.size());
                            static_assert(sizeof(opss[tid].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
                            static_assert(sizeof(fss[tid].ids()[0]) == sizeof(uint64_t), "Must be 64-bit");
                            ptr = opss.size() ? opss[tid].ids().data(): fss[tid].ids().data();
                            checked_fwrite(ofp, ptr, sizeof(uint64_t) * ss);
                        } else {
                            std::fprintf(stderr, "sspace: %d. result type: %d\n", opts.sspace_, opts.kmer_result_);
                            throw std::runtime_error("Not PSET, MULTISET, ONEPERM/FULLSETSKETCH.");
                        }
                        if(ptr && opts.build_mmer_matrix_) {
                            checked_fwrite(ofp, ptr, sizeof(uint64_t) * ss);
                            std::fprintf(stderr, "About to copy to kmers of size %zu\n", ret.kmers_.size());
                            std::copy(ptr, ptr + ss, &ret.kmers_[i * ss]);
                        } else std::fprintf(stderr, "mmer matrix not built\n");
                        if(ofp) std::fclose(ofp);
                    }
                    //std::fprintf(stderr, "Save kmercounts\n");
                    if(opts.save_kmercounts_ || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
                        assert(ret.kmercountfiles_.size());
                        DBG_ONLY(std::fprintf(stderr, "About to save kmer counts to %s\n", destkmercounts.data());)
                        ret.kmercountfiles_[i] = destkmercounts;
                        if((ofp = std::fopen(destkmercounts.data(), "wb")) == nullptr) throw std::runtime_error("Failed to write k-mer counts");
                        std::vector<float> tmp(ss);
                        if(opts.kmer_result_ == FULL_MMER_COUNTDICT || (opts.kmer_result_ == FULL_MMER_SET && opts.save_kmercounts_)) {
                            tmp.resize(kmerveccounts.size());
                            std::copy(kmerveccounts.data(), kmerveccounts.data() + kmerveccounts.size(), tmp.data());
                        } else if(pmhs.size()) {
                            auto p = pmhs[tid].idcounts().data();
                            std::copy(p, p + pmhs[tid].idcounts().size(), tmp.data());
                        } else if(bmhs.size()) {
                            auto p = bmhs[tid].idcounts().data();
                            std::copy(p, p + bmhs[tid].idcounts().size(), tmp.data());
                        } else if(opss.size()) {
                            auto p = opss[tid].idcounts().data();
                            std::copy(p, p + opss[tid].idcounts().size(), tmp.data());
                        } else if(fss.size()) {
                            auto p = fss[tid].idcounts().data();
                            std::copy(p, p + fss[tid].idcounts().size(), tmp.data());
                        }
                        std::fprintf(stderr, "Copying from tmp to fp (%p)\n", (void *)ofp);
                        checked_fwrite(ofp, tmp.data(), tmp.size() * sizeof(float));
                        tmp.resize(ss);
                        std::fclose(ofp);
                        if(ret.kmercounts_.size()) {
                            std::fprintf(stderr, "Copying range of size %zu from tmp to ret.kmercounts of size %zu\n", tmp.size(), ret.kmercounts_.size());
                            std::copy(tmp.begin(), tmp.end(), &ret.kmercounts_[i * ss]);
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
                            std::fprintf(stderr, "Pushing back. l%zu, newm%zu\n", l, newm);
                            void *newptr = nullptr;
                            if(posix_memalign((void **)&newptr, 16, newm * sizeof(DT))) throw std::bad_alloc();
                            assert(newptr);
                            std::copy((DT *)dptr, (DT *)dptr + m, (DT *)newptr);
                            dptr = newptr;ptr = (DT *)dptr;
                            m = newm;
                        }
                        if(opts.homopolymer_compress_minimizers_) {
                            std::fprintf(stderr, "hpc\n");
                            if(l > 0 && ptr[l - 1] == x) return;
                        }
                        ptr[l++] = x;
                    });
                    assert(dptr);
                    checked_fwrite(ofp, dptr, l * (1 + opts.use128()) * sizeof(uint64_t));
                    std::free(dptr);
                    std::fclose(ofp);
                } else if(opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) {
                    // These occur twice, because if the user asks for counts, or if the user asks for a minimum count level for inclusion.
                    // Because of this, we have to generate the key-count map.
                    // Those cases are handled above with the count-based methods.
                    std::fprintf(stderr, "kmer result is oneperm or setsketch\n");
                    std::FILE * ofp;
                    if((ofp = std::fopen(destination.data(), "wb")) == nullptr) throw std::runtime_error("Failed to open file for writing minimizer sequence");
                    if(opss.empty() && fss.empty()) throw std::runtime_error("Both opss and fss are empty\n");
                    const size_t opsssz = opss.size();
                    perf_for_substrs([&](auto hv) {
                        if(opsssz) opss.at(tid).update(hv);
                        else fss.at(tid).update(hv);
                    });
                    const uint64_t *ids = nullptr;
                    const uint32_t *counts = nullptr;
                    const RegT *ptr = opsssz ? opss[tid].data(): fss[tid].data();
                    assert(ptr);
                    if(opts.build_mmer_matrix_)
                        ids = opsssz ? opss[tid].ids().data(): fss[tid].ids().data();
                    if(opts.build_count_matrix_) {
                        counts = opsssz ? opss[tid].idcounts().data(): fss[tid].idcounts().data();
                    }
                    checked_fwrite(ofp, ptr, sizeof(RegT) * ss);
                    std::fclose(ofp);
                    if(ptr && ret.signatures_.size()) std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                    if(ids && ret.kmers_.size())
                        std::copy(ids, ids + ss, &ret.kmers_[i * ss ]);
                    if(counts && ret.kmercounts_.size())
                        std::copy(counts, counts + ss, &ret.kmercounts_[i * ss]);
                } else throw std::runtime_error("Unexpected: Not FULL_MMER_SEQUENCE, FULL_MMER_SET, ONE_PERM, FULL_SETSKETCH, SPACE_MULTISET, or SPACE_PSET");
            }
        }
    }
    for(size_t i = 0; i < nt; ++i) bns::kseq_destroy_stack(kseqs[i]);
    std::free(kseqs);
    return ret;
}

} // dashing2
