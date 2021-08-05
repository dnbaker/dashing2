#include "fastxsketch.h"

namespace dashing2 {
using bns::InputType;


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
    bns::InputType it_;
    OptSketcher&operator=(const OptSketcher &o) {
        rh_ = o.rh_; rh128_ = o.rh128_; enc_ = o.enc_; enc128_ = o.enc128_; k_ = o.k_; w_ = o.w_; use128_ = o.use128_; it_ = o.it_;
        input_mode(o.input_mode());
        if(o.ctr) ctr.reset(new Counter(*o.ctr));
        if(o.bmh) bmh.reset(new BagMinHash(*o.bmh));
        else if(o.pmh) pmh.reset(new ProbMinHash(*o.pmh));
        else if(o.opss) opss.reset(new OPSetSketch(o.opss->size()));
        else if(o.fss) fss.reset(new FullSetSketch(*o.fss));
        else if(o.omh) omh.reset(new OrderMinHash(*o.omh));
        return *this;
    }
    OptSketcher(const OptSketcher &o): rh_(o.rh_), rh128_(o.rh128_), enc_(o.enc_), enc128_(o.enc128_), k_(o.k_), w_(o.w_), use128_(o.use128_) {
        input_mode(o.input_mode());
        if(o.ctr) ctr.reset(new Counter(*o.ctr));
        if(o.bmh) bmh.reset(new BagMinHash(*o.bmh));
        else if(o.pmh) pmh.reset(new ProbMinHash(*o.pmh));
        else if(o.opss) opss.reset(new OPSetSketch(o.opss->size()));
        else if(o.fss) fss.reset(new FullSetSketch(*o.fss));
        else if(o.omh) omh.reset(new OrderMinHash(*o.omh));
    }
    OptSketcher(const Dashing2Options &opts): enc_(opts.enc_), enc128_(std::move(enc_.to_u128())), k_(opts.k_), w_(opts.w_), use128_(opts.use128()) {
        input_mode(opts.input_mode());
        assert(opts.hashtype() == it_);
        if(opts.sspace_ == SPACE_EDIT_DISTANCE)
            omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
    }
    template<typename Func>
    void for_each(const Func &func, const char *s, size_t n) {
        if(use128_) {
            if(unsigned(k_) <= enc_.nremperres128()) {
                if(entmin)
                    enc128_.to_entmin128().for_each(func, s, n);
                else
                    enc128_.for_each(func, s, n);
            } else {
                rh128_.for_each(func, s, n);
            }
        } else {
            if(unsigned(k_) <= enc_.nremperres64()) {
                if(entmin)
                    enc_.to_entmin64().for_each(func, s, n);
                else
                    enc_.for_each(func, s, n);
            } else
                rh_.for_each_hash(func, s, n);
        }
    }
    void input_mode(InputType it) {it_ = it; enc_.hashtype(it); enc128_.hashtype(it); rh_.hashtype(it); rh128_.hashtype(it);}
    bns::InputType input_mode() const {return it_;}
    bool enable_protein() const {return it_;}
    void reset() {
        rh_.reset(); rh128_.reset();
        if(ctr) ctr->reset();
        if(bmh) bmh->reset();
        else if(pmh) pmh->reset();
        else if(opss) opss->reset();
        else if(fss) fss->reset();
        //if(omh) omh->reset();
    }
};
void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, std::vector<OptSketcher> &sketchvec, size_t &lastindex, size_t nthreads);

FastxSketchingResult fastx2sketch_byseq(Dashing2Options &opts, const std::string &path, kseq_t *kseqs, bool parallel, size_t seqs_per_batch) {
    const bool save_ids = opts.save_kmers_ || opts.build_mmer_matrix_;
    const bool save_idcounts = opts.save_kmercounts_ || opts.build_count_matrix_;
    FastxSketchingResult ret;
    gzFile ifp;
    kseq_t *myseq = kseqs ? &kseqs[OMP_ELSE(omp_get_thread_num(), 0)]: (kseq_t *)std::calloc(sizeof(kseq_t), 1);
    size_t batch_index = 0;
    OptSketcher sketcher(opts);
    sketcher.rh_ = opts.rh_;
    sketcher.rh128_ = opts.rh128_;
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm()) {
            //std::fprintf(stderr, "Using OPSS\n");
            sketcher.opss.reset(new OPSetSketch(opts.sketchsize_));
            if(opts.count_threshold_ > 0) sketcher.opss->set_mincount(opts.count_threshold_);
        } else sketcher.fss.reset(new FullSetSketch(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_MULTISET) {
            //std::fprintf(stderr, "Using BMH\n");
        sketcher.bmh.reset(new BagMinHash(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_PSET) {
        sketcher.pmh.reset(new ProbMinHash(opts.sketchsize_));
        //std::fprintf(stderr, "Setting sketcher.pmh: %p\n", (void *)sketcher.pmh.get());
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        sketcher.omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
        //std::fprintf(stderr, "Setting sketcher.omh: %p\n", (void *)sketcher.omh.get());
    } else THROW_EXCEPTION(std::runtime_error("Should have been set space, multiset, probset, or edit distance"));
    if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET) {
        sketcher.ctr.reset(new Counter(opts.cssize()));
    }

    std::vector<OptSketcher> sketching_data;
    size_t nt = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        nt = omp_get_num_threads();
    }
#endif
    if(!parallel) nt = 1;
    std::vector<OptSketcher> sv;
    if(nt > 1) {
        sketching_data.resize(nt, sketcher);
    } else sketching_data.emplace_back(std::move(sketcher));
    //std::fprintf(stderr, "save ids: %d, save counts %d\n", save_ids, save_idcounts);
    size_t lastindex = 0;
    for_each_substr([&](const auto &x) {
        if((ifp = gzopen(x.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to read from ") + x));
        gzbuffer(ifp, 1u << 17);
        bns::kseq_assign(myseq, ifp);
        for(int c;(c = kseq_read(myseq)) >= 0;) {
            //std::fprintf(stderr, "Sequence %s of length %zu\n", myseq->name.s, myseq->seq.l);
            ret.sequences_.emplace_back(myseq->seq.s, myseq->seq.l);
            const int off = (myseq->name.s[0] == '>');
            ret.names_.emplace_back(myseq->name.s + off, myseq->name.l - off);
            if(++batch_index == seqs_per_batch) {
                //std::fprintf(stderr, "batch index = %zu\n", batch_index);
                resize_fill(opts, ret, seqs_per_batch, sketching_data, lastindex, nt);
                batch_index = 0;
                seqs_per_batch = std::min(seqs_per_batch << 1, size_t(0x1000));
            }
        }
        gzclose(ifp);
    }, path);
    if(!kseqs) kseq_destroy(myseq);
    if(batch_index) resize_fill(opts, ret, batch_index, sketching_data, lastindex, nt);
    ret.names_.resize(lastindex);
    ret.sequences_.resize(lastindex);
    ret.cardinalities_.resize(lastindex);
    if(ret.kmers_.size()) ret.kmers_.resize(opts.sketchsize_ * lastindex);
    if(ret.kmercounts_.size()) ret.kmercounts_.resize(opts.sketchsize_ * lastindex);
#ifndef NDEBUG
    std::fprintf(stderr, "ret kmer size %zu\n", ret.kmers_.size());
    std::fprintf(stderr, "ret names size %zu\n", ret.names_.size());
#endif
    return ret;
}
void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, std::vector<OptSketcher> &sketchvec, size_t &lastindex, size_t nt) {
    //std::fprintf(stderr, "Calling resize_fill with newsz = %zu\n", newsz);
    const size_t oldsz = ret.names_.size();
    newsz = oldsz + newsz;
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE && opts.build_sig_matrix_) {
        //std::fprintf(stderr, "old sig size %zu, new %zu\n", ret.signatures_.size(), newsz * opts.sketchsize_);
        ret.signatures_.resize(newsz * opts.sketchsize_);
    }
    ret.cardinalities_.resize(newsz);
    //std::fprintf(stderr, "mmer matrix size %zu. buuild %d, save kmers %d\n", ret.kmers_.size(), opts.build_mmer_matrix_, opts.save_kmers_);
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE && (opts.build_mmer_matrix_ || opts.save_kmers_)) {
        ret.kmers_.resize(opts.sketchsize_ * newsz);
    }
    if((opts.kmer_result_ != FULL_MMER_SEQUENCE) && (opts.build_count_matrix_ || opts.save_kmercounts_)) {
        ret.kmercounts_.resize(opts.sketchsize_ * newsz);
    }
    //std::fprintf(stderr, "Parsing %s\n", sketchvec.front().enable_protein() ? "Protein": "DNA");
    std::unique_ptr<std::vector<uint64_t>[]> seqmins;
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
        seqmins.reset(new std::vector<uint64_t>[(oldsz - lastindex)]);
    }
    OMP_PRAGMA("omp parallel for num_threads(nt) schedule(dynamic,16)")
    for(size_t i = lastindex; i < oldsz; ++i) {
        const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        DBG_ONLY(std::fprintf(stderr, "%zu/%zu -- parsing sequence from tid = %d\n", i, oldsz, tid);)
        auto &sketchers(sketchvec[tid]);
        sketchers.reset();
        if(sketchers.omh) {
            std::vector<uint64_t> res = sketchers.omh->hash(ret.sequences_[i].data(), ret.sequences_[i].size());
            auto destp = &ret.signatures_[i * opts.sketchsize_];
            if constexpr(sizeof(RegT) == sizeof(uint64_t)) { // double
                // Copy as-is
                std::memcpy(destp, res.data(), res.size() * sizeof(uint64_t));
            } else if constexpr(sizeof(RegT) > sizeof(uint64_t)) { // long double
                std::transform(res.begin(), res.end(), (u128_t *)destp, [](uint64_t x) {
                    long double ret;
                    std::memcpy(&ret, &x, sizeof(x));
                    ((uint64_t *)&ret)[1] = wy::wyhash64_stateless(&x);
                    return ret;
                });
            } else {
                std::copy(res.begin(), res.end(), (uint32_t *)destp);
            }
            ret.cardinalities_[i] = ret.sequences_[i].size();
        } else if(seqmins) {
            auto &myseq(seqmins[i - lastindex]);
            sketchers.for_each([&](auto x) {
                x = maskfn(x);
                if(opts.fs_ && opts.fs_->in_set(x)) return;
                if(!opts.homopolymer_compress_minimizers_ || myseq.empty() ||
                    (sizeof(x) == 8 ? myseq.back() != x: std::memcmp(&myseq[myseq.size() - 2], &x, 16) != 0))
                {
                    if(sizeof(x) == 8) myseq.push_back(x);
                    else {
                        const size_t oldsz = myseq.size();
                        myseq.resize(oldsz + 2);
                        std::memcpy(&myseq[oldsz], &x, 16);
                    }
                }
            }, ret.sequences_[i].data(), ret.sequences_[i].size());
            // This handles the case where the sequence is shorter than the window size
            // and the entire sequence yields no minimizer.
            // We instead take the minimum-hashed value from the b-tree
            if(opts.w_ > opts.k_ && myseq.empty() && ret.sequences_[i].size()) {
                //std::fprintf(stderr, "Adding in single item for small seq]\n");
                if(sketchers.rh128_.n_in_queue()) {
                    u128_t v = sketchers.rh128_.max_in_queue().el_;
                    myseq.push_back(0); myseq.push_back(0);
                    std::memcpy(&myseq[myseq.size() - 2], &v, sizeof(v));
                } else if(sketchers.enc128_.n_in_queue()) {
                    u128_t v = sketchers.enc128_.max_in_queue().el_;
                    myseq.push_back(0); myseq.push_back(0);
                    std::memcpy(&myseq[myseq.size() - 2], &v, sizeof(v));
                } else if(sketchers.rh_.n_in_queue()) {
                    myseq.push_back(sketchers.rh_.max_in_queue().el_);
                } else if(sketchers.enc_.n_in_queue()) {
                    myseq.push_back(sketchers.enc_.max_in_queue().el_);
                }
            }
            ret.cardinalities_[i] = myseq.size();
            //std::fprintf(stderr, "Processed items for %zu\n", i);
        } else {
            //std::fprintf(stderr, "Calcing hash for seq = %zu/%s\n", i,  ret.sequences_[i].data());
            auto seqp = ret.sequences_[i].data();
            auto seql = ret.sequences_[i].size();
            if(opts.fs_) {
                sketchers.for_each([&](auto x) {
                    x = maskfn(x);
                    if(opts.fs_->in_set(x)) return;
                    if(sketchers.opss) sketchers.opss->update(x);
                    else if(sketchers.ctr)
                        sketchers.ctr->add(x);
                    else if(sketchers.fss)
                        sketchers.fss->update(x);
                }, seqp, seql);
            } else {
                const bool isop = sketchers.opss.get(), isctr = sketchers.ctr.get(), isfs = sketchers.fss.get();
                sketchers.for_each([&](auto x) {
                    x = maskfn(x);
                    if(isop) sketchers.opss->update(x);
                    else if(isctr) sketchers.ctr->add(x);
                    else if(isfs) sketchers.fss->update(x);
                }, seqp, seql);
            }
            RegT *ptr = nullptr;
            const uint64_t *kmer_ptr = nullptr;
            std::vector<double> kmercounts;
            if(opts.sspace_ == SPACE_SET) {
                if(sketchers.ctr) {
                    if(sketchers.opss) {
                        sketchers.ctr->finalize(*sketchers.opss, opts.count_threshold_);
                        ptr = sketchers.opss->data();
                        ret.cardinalities_[i] = sketchers.opss->getcard();
                    } else {
                        sketchers.ctr->finalize(*sketchers.fss, opts.count_threshold_);
                        ptr = sketchers.fss->data();
                        ret.cardinalities_[i] = sketchers.fss->getcard();
                    }
                } else {
                    ptr = sketchers.opss ? sketchers.opss->data(): sketchers.fss->data();
                    ret.cardinalities_[i] = sketchers.opss ? sketchers.opss->getcard(): sketchers.fss->getcard();
                    if(ret.cardinalities_[i] < 10 * opts.sketchsize_) {
                        flat_hash_set<uint64_t> ids;
                        ids.reserve(opts.sketchsize_);
                        if(opts.fs_) {
                            sketchers.for_each([&](auto x) {
                                x = maskfn(x);
                                if(opts.fs_->in_set(x)) return;
                                ids.insert(x);
                            }, ret.sequences_[i].data(), ret.sequences_[i].size());
                        } else {
                            sketchers.for_each([&](auto x) {ids.insert(maskfn(x));}, ret.sequences_[i].data(), ret.sequences_[i].size());
                        }
                        ret.cardinalities_[i] = ids.size();
                    }
                    kmer_ptr = sketchers.opss ? sketchers.opss->ids().data(): sketchers.fss->ids().data();
                    kmercounts.resize(opts.sketchsize_);
                    if(sketchers.opss) {
                        auto &idc = sketchers.opss->idcounts();
                        std::copy(idc.begin(), idc.end(), kmercounts.begin());
                    } else if(sketchers.fss->idcounts().size()) {
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
            } else THROW_EXCEPTION(std::runtime_error("Not yet implemented?"));
            std::copy(ptr, ptr + opts.sketchsize_, &ret.signatures_[i * opts.sketchsize_]);
            if(kmer_ptr && ret.kmers_.size()) {
                //std::fprintf(stderr, "Copying k-mers out, kmers size %zu, idx = %zu\n", ret.kmers_.size(), i * opts.sketchsize_);
                std::copy(kmer_ptr, kmer_ptr + opts.sketchsize_, &ret.kmers_[i * opts.sketchsize_]);
            }// else std::fprintf(stderr, "k-mers not saved\n");
            if(kmercounts.size() && ret.kmercounts_.size()) {
                std::copy(kmercounts.begin(), kmercounts.end(), &ret.kmercounts_[i * opts.sketchsize_]);
            }
        }
    }
    if(seqmins) {
        const size_t seqminsz = oldsz - lastindex;
        using OT = std::conditional_t<(sizeof(RegT) == 8), uint64_t, std::conditional_t<(sizeof(RegT) == 4), uint32_t, u128_t>>;
        static_assert(sizeof(RegT) == sizeof(OT), "Size of hash registers must match that being sketched");
        for(size_t i = 0; i < seqminsz; ++i) {
            //std::fprintf(stderr, "Copying out items from %zu/%zu\n", i, seqminsz);
            const auto &x(seqmins[i]);
            const OT *ptr = (const OT *)x.data();
            size_t xsz = x.size();
            if constexpr(sizeof(RegT) == 16) xsz >>= 1;
            ret.nperfile_.push_back(xsz);
            if constexpr(sizeof(RegT) == 4)
                std::copy((const uint64_t *)ptr, (const uint64_t *)ptr + xsz, std::back_inserter(ret.signatures_));
            else
                ret.signatures_.insert(ret.signatures_.end(), ptr, ptr + xsz);
        }
    }
    lastindex = oldsz;
}
} // dashing2
