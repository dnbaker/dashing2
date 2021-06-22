#include "fastxsketch.h"

namespace dashing2 {

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
    bool use128_, enable_protein_;
    OptSketcher(const Dashing2Options &opts): enc_(opts.enc_), enc128_(std::move(enc_.to_u128())), k_(opts.k_), w_(opts.w_), use128_(opts.use128()), enable_protein_(opts.parse_protein()) {
        if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
            omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
        }
    }
    template<typename Func>
    void for_each(const Func &func, const char *s, size_t n) {
        if(k_ > 64 || enable_protein_) {
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
    bool enable_protein() const {return enable_protein_;}
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
void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, OptSketcher &sketchers, size_t &lastindex);

FastxSketchingResult fastx2sketch_byseq(Dashing2Options &opts, const std::string &path, kseq_t *kseqs, size_t seqs_per_batch) {
    FastxSketchingResult ret;
    gzFile ifp;
    kseq_t *myseq = kseqs ? &kseqs[OMP_ELSE(omp_get_thread_num(), 0)]: (kseq_t *)std::calloc(sizeof(kseq_t), 1);
    size_t batch_index = 0;
    OptSketcher sketching_data(opts);
    sketching_data.rh_ = opts.rh_;
    sketching_data.rh128_ = opts.rh128_;
    const bool save_ids = opts.save_kmers_ || opts.build_mmer_matrix_;
    const bool save_idcounts = opts.save_kmercounts_ || opts.build_count_matrix_;
    //std::fprintf(stderr, "save ids: %d, save counts %d\n", save_ids, save_idcounts);
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm()) {
            //std::fprintf(stderr, "Using OPSS\n");
            sketching_data.opss.reset(new OPSetSketch(opts.sketchsize_));
            if(opts.count_threshold_ > 0) sketching_data.opss->set_mincount(opts.count_threshold_);
        } else sketching_data.fss.reset(new FullSetSketch(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_MULTISET) {
            std::fprintf(stderr, "Using BMH\n");
        sketching_data.bmh.reset(new BagMinHash(opts.sketchsize_, save_ids, save_idcounts));
    } else if(opts.sspace_ == SPACE_PSET) {
        sketching_data.pmh.reset(new ProbMinHash(opts.sketchsize_));
        std::fprintf(stderr, "Setting sketching_data.pmh: %p\n", (void *)sketching_data.pmh.get());
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        sketching_data.omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
        std::fprintf(stderr, "Setting sketching_data.omh: %p\n", (void *)sketching_data.omh.get());
    } else throw std::runtime_error("Should have been set space, multiset, probset, or edit distance");
    if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET || (opts.sspace_ == SPACE_SET && opts.kmer_result_ == FULL_SETSKETCH && opts.count_threshold_ > 0.)) {
        sketching_data.ctr.reset(new Counter(opts.cssize()));
    }

    size_t lastindex = 0;
    for_each_substr([&](const auto &x) {
        if((ifp = gzopen(x.data(), "rb")) == nullptr) throw std::runtime_error(std::string("Failed to read from ") + x);
        gzbuffer(ifp, 1u << 17);
        bns::kseq_assign(myseq, ifp);
        for(int c;(c = kseq_read(myseq)) >= 0;) {
            std::fprintf(stderr, "Sequence %s of length %zu\n", myseq->name.s, myseq->seq.l);
            ret.sequences_.emplace_back(myseq->seq.s, myseq->seq.l);
            const int off = (myseq->name.s[0] == '>');
            ret.names_.emplace_back(myseq->name.s + off, myseq->name.l - off);
            if(++batch_index == seqs_per_batch) {
                std::fprintf(stderr, "batch index = %zu\n", batch_index);
                resize_fill(opts, ret, seqs_per_batch, sketching_data, lastindex);
                batch_index = 0;
                seqs_per_batch = std::min(seqs_per_batch << 1, size_t(0x1000));
            }
        }
        gzclose(ifp);
    }, path);
    if(!kseqs) kseq_destroy(myseq);
    if(batch_index) resize_fill(opts, ret, batch_index, sketching_data, lastindex);
    ret.names_.resize(lastindex);
    ret.sequences_.resize(lastindex);
    ret.cardinalities_.resize(lastindex);
    if(ret.kmers_.size()) ret.kmers_.resize(opts.sketchsize_ * lastindex);
    if(ret.kmercounts_.size()) ret.kmercounts_.resize(opts.sketchsize_ * lastindex);
    std::fprintf(stderr, "ret kmer size %zu\n", ret.kmers_.size());
    std::fprintf(stderr, "ret names size %zu\n", ret.names_.size());
    return ret;
}
void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, OptSketcher &sketchers, size_t &lastindex) {
    std::fprintf(stderr, "Calling resize_fill with newsz = %zu\n", newsz);
    const size_t oldsz = ret.names_.size();
    newsz = oldsz + newsz;
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE && opts.build_sig_matrix_) {
        //std::fprintf(stderr, "old sig size %zu, new %zu\n", ret.signatures_.size(), newsz * opts.sketchsize_);
        ret.signatures_.resize(newsz * opts.sketchsize_);
    }
    ret.cardinalities_.resize(newsz);
    //std::fprintf(stderr, "mmer matrix size %zu. buuild %d, save kmers %d\n", ret.kmers_.size(), opts.build_mmer_matrix_, opts.save_kmers_);
    if(opts.build_mmer_matrix_ || opts.save_kmers_) {
        ret.kmers_.resize(opts.sketchsize_ * newsz);
    }
    if(opts.build_count_matrix_ || opts.save_kmercounts_) {
        ret.kmercounts_.resize(opts.sketchsize_ * newsz);
    }
    std::fprintf(stderr, "Parsing %s\n", sketchers.enable_protein() ? "Protein": "DNA");
    std::fprintf(stderr, "About to go through list, %zu sigsize, %zu mmer size %zu mmer count size\n", ret.signatures_.size(), ret.kmers_.size(), ret.kmercounts_.size());
    std::unique_ptr<std::vector<uint64_t>[]> seqmins;
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
        seqmins.reset(new std::vector<uint64_t>[(oldsz - lastindex)]);
    }
    for(size_t i = lastindex; i < oldsz; ++i) {
        std::fprintf(stderr, "%zu/%zu -- hashing sequence\n", i, newsz);
        //std::fprintf(stderr, "Seq %zu/%s\n", seqp.second, seqp.first);
        sketchers.reset();
        if(sketchers.omh) {
            std::fprintf(stderr, "omh is set, not hashing: %s/%zu\n", ret.sequences_[i].data(), ret.sequences_[i].size());
            std::vector<uint64_t> res = sketchers.omh->hash(ret.sequences_[i].data(), ret.sequences_[i].size());
            auto destp = &ret.signatures_[i * opts.sketchsize_];
            if constexpr(sizeof(RegT) == sizeof(uint64_t)) { // double
                std::fprintf(stderr, "Copying 64-bit registers for OMH\n");
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
                if(opts.fs_ && opts.fs_->in_set(x)) return;
                if(myseq.empty() || !opts.homopolymer_compress_minimizers_ || myseq.back() != x)
                    myseq.emplace_back(x);
            }, ret.sequences_[i].data(), ret.sequences_[i].size());
        } else {
            assert(!sketchers.opss || sketchers.opss->total_updates() == 0u);
            std::fprintf(stderr, "Calcing hash for seq = %zu/%s\n", i,  ret.sequences_[i].data());
            sketchers.for_each([&](auto x) {
                if(opts.fs_ && opts.fs_->in_set(x)) return;
                if(sketchers.opss) sketchers.opss->update(x);
                else if(sketchers.ctr)
                    sketchers.ctr->add(x);
                else if(sketchers.fss)
                    sketchers.fss->update(x);
            }, ret.sequences_[i].data(), ret.sequences_[i].size());
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
                } else if(sketchers.opss) {
                    ptr = sketchers.opss->data();
                    ret.cardinalities_[i] = sketchers.opss->getcard();
                    if(ret.cardinalities_[i] < 10 * opts.sketchsize_) {
                        ska::flat_hash_set<uint64_t> ids;
                        ids.reserve(opts.sketchsize_);
                        sketchers.for_each([&](auto x) {ids.insert(x);}, ret.sequences_[i].data(), ret.sequences_[i].size());
                        ret.cardinalities_[i] = ids.size();
                    }
                    kmer_ptr = sketchers.opss->ids().data();
                    kmercounts.resize(opts.sketchsize_);
                    auto &idc = sketchers.opss->idcounts();
                    std::copy(idc.begin(), idc.end(), kmercounts.begin());
                } else {
                    assert(sketchers.fss);
                    ptr = sketchers.fss->data();
                    ret.cardinalities_[i] = sketchers.fss->getcard();
                    DBG_ONLY(std::fprintf(stderr, "Card is %g from FSS\n", sketchers.fss->getcard());)
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
            std::fprintf(stderr, "Copying from %p to container of size %zu at %zu\n", (void *)ptr, ret.signatures_.size(), size_t(i * opts.sketchsize_));
            std::copy(ptr, ptr + opts.sketchsize_, &ret.signatures_[i * opts.sketchsize_]);
            if(kmer_ptr && ret.kmers_.size()) {
                std::fprintf(stderr, "Copying k-mers out, kmers size %zu, idx = %zu\n", ret.kmers_.size(), i * opts.sketchsize_);
                std::copy(kmer_ptr, kmer_ptr + opts.sketchsize_, &ret.kmers_[i * opts.sketchsize_]);
            } else std::fprintf(stderr, "k-mers not saved\n");
            if(kmercounts.size() && ret.kmercounts_.size()) {
                std::fprintf(stderr, "Copying k-mer counts out\n");
                std::copy(kmercounts.begin(), kmercounts.end(), &ret.kmercounts_[i * opts.sketchsize_]);
                std::fprintf(stderr, "Copied k-mer counts out\n");
            }
        }
    }
    if(seqmins) {
        const size_t seqminsz = oldsz - lastindex;
        static_assert(sizeof(RegT) == 8, "Must use doubles, sorry");
        for(size_t i = 0; i < seqminsz; ++i) {
            const auto &x(seqmins[i]);
            const uint64_t *ptr = (const uint64_t *)x.data();
            ret.nperfile_.push_back(x.size());
            ret.signatures_.insert(ret.signatures_.end(), ptr, ptr + x.size());
        }
    }
    lastindex = oldsz;
}
} // dashing2
