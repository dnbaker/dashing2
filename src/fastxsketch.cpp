#include "d2.h"
#include <optional>
namespace dashing2 {

struct FastxSketchingResult {
    std::vector<std::string> names_;
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
        if(opts.one_perm_) make_save(opss);
        else make_save(fss);
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
    if(!opts.parse_by_seq_) {
        ret.nperfile_.resize(paths.size());
        std::fprintf(stderr, "Parse by seq not available\n");
        std::exit(1);
    } else {
        if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            throw std::runtime_error("edit distance is only available in parse by seq mode");
        }
        if(paths.size() == 1) {
            // Parallelized construction
            throw std::runtime_error("Not yet completed: parallel sketching of a single file");
        } else {
            // TODO:
            // Change ordering to process larger datasets first for better load-balancing
            ret.names_.resize(paths.size());
            std::copy(paths.begin(), paths.end(), ret.names_.begin());
            ret.signatures_.resize(ss * paths.size());
            if(opts.save_kmers_) {
                ret.kmers_.resize(ss * paths.size());
                if(opts.save_kmercounts_)
                    ret.kmercounts_.resize(ss * paths.size());
            }
            OMP_PFOR_DYN
            for(size_t i = 0; i < paths.size(); ++i) {
                const int tid = OMP_ELSE(omp_get_thread_num(), 0);
                kseq_t *ksp = kseqs + tid;
                auto &path = paths[i];
                reset(tid);
                auto perf_for_substrs = [&](const auto &func) {
                    for_each_substr([&](const std::string &subpath) {
                        if(opts.use128()) {
                            opts.rh128_.for_each_hash(func, subpath.data(), ksp);
                        } else {
                            opts.rh_.for_each_hash(func, subpath.data(), ksp);
                        }
                    }, path);
                };
                switch(opts.sspace_) {
                    case SPACE_SET: {
                        const RegT *ptr;
                        const uint64_t *ids;
                        const uint32_t *counts;
#define SPACE_SET_BODY\
                            perf_for_substrs([&](auto hv) {s.update(hv);});\
                            ptr = s.data();\
                            ids = s.ids().data();\
                            if(opts.save_kmercounts_) counts = s.idcounts().data();
                        if(opts.one_perm_) {
                            auto &s = opss[tid];
                            SPACE_SET_BODY
                        } else {
                            auto &s = fss[tid];
                            SPACE_SET_BODY
#undef SPACE_SET_BODY
                        }
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                        if(opts.save_kmers_) {
                            std::copy(ids, ids + ss, &ret.kmers_[i * ss ]);
                            if(opts.save_kmercounts_) {
                                std::copy(counts, counts + ss, &ret.kmercounts_[i * ss]);
                            }
                        }
                    }
                    break;
                    case SPACE_MULTISET: case SPACE_PSET: {
                        // We are doing summation currently, and sketching will be performed afterward
                        auto &ctr = ctrs[tid];
                        perf_for_substrs([&](auto x) {ctr.add(x);});
                        RegT *ptr;
                        if(opts.sspace_ == SPACE_MULTISET) {
                            ctr.finalize(bmhs[tid]);
                            ptr = bmhs[tid].data();
                        } else {
                            ctr.finalize(pmhs[tid]);
                            ptr = pmhs[tid].data();
                        }
                        std::copy(ptr, ptr + ss, &ret.signatures_[i * ss]);
                        if(opts.save_kmers_) {
                            if(opts.sspace_ == SPACE_MULTISET) {
                                std::fprintf(stderr, "Not yet supported: k-mer id saving from BagMinHash (multiset sketching).\n");
                            } else {
                                std::copy(pmhs[tid].res_.begin(), pmhs[tid].res_.end(), &ret.kmers_[i * ss]);
                                if(opts.save_kmercounts_) {
                                    std::copy(pmhs[tid].rcounts_.begin(), pmhs[tid].rcounts_.end(), &ret.kmercounts_[i * ss]);
                                }
                            }
                        }
                    }
                    break;
                    default: __builtin_unreachable();
                        throw std::invalid_argument("Unexpected space: not set, multiset, pset, or edit distance");
                }
            }
        }
    }
    for(size_t i = 0; i < nt; ++i) bns::kseq_destroy_stack(kseqs[i]);
    std::free(kseqs);
    return ret;
}

} // dashing2
