#include "fastxsketch.h"

namespace dashing2 {
using namespace variation;
using bns::InputType;

#undef DBG_ONLY
#define DBG_ONLY(x) do {if(verbosity >= DEBUG) {x}} while(0)

struct OptSketcher {
    std::unique_ptr<BagMinHash> bmh;
    std::unique_ptr<ProbMinHash> pmh;
    std::unique_ptr<OPSetSketch> opss;
    std::unique_ptr<FullSetSketch> fss;
    std::unique_ptr<OrderMinHash> omh;
    std::unique_ptr<VSetSketch> cfss;
    std::unique_ptr<Counter> ctr;
    bns::RollingHasher<uint64_t> rh_;
    bns::RollingHasher<u128_t> rh128_;
    bns::Encoder<bns::score::Lex, uint64_t> enc_;
    bns::Encoder<bns::score::Lex, u128_t> enc128_;
    bns::Encoder<bns::score::Entropy, uint64_t> ence_;
    bns::Encoder<bns::score::Entropy, u128_t> ence128_;
    int k_, w_;
    bool use128_;
    bns::InputType it_;
    OptSketcher&operator=(const OptSketcher &o) {
        rh_ = o.rh_; rh128_ = o.rh128_; enc_ = o.enc_; enc128_ = o.enc128_; ence_ = o.ence_; ence128_ = o.ence128_; k_ = o.k_; w_ = o.w_; use128_ = o.use128_; it_ = o.it_;
        input_mode(o.input_mode());
        if(o.ctr) ctr.reset(new Counter(*o.ctr));
        if(o.bmh) bmh.reset(new BagMinHash(*o.bmh));
        else if(o.pmh) pmh.reset(new ProbMinHash(*o.pmh));
        else if(o.opss) opss.reset(new OPSetSketch(o.opss->size()));
        else if(o.fss) fss.reset(new FullSetSketch(*o.fss));
        else if(o.omh) omh.reset(new OrderMinHash(*o.omh));
        else if(o.cfss) cfss.reset(new VSetSketch(*o.cfss));
        return *this;
    }
    OptSketcher(const OptSketcher &o): rh_(o.rh_), rh128_(o.rh128_), enc_(o.enc_), enc128_(o.enc128_), ence_(o.ence_), ence128_(o.ence128_), k_(o.k_), w_(o.w_), use128_(o.use128_) {
        input_mode(o.input_mode());
        if(o.ctr) ctr.reset(new Counter(*o.ctr));
        if(o.bmh) bmh.reset(new BagMinHash(*o.bmh));
        else if(o.pmh) pmh.reset(new ProbMinHash(*o.pmh));
        else if(o.opss) opss.reset(new OPSetSketch(o.opss->size()));
        else if(o.fss) fss.reset(new FullSetSketch(*o.fss));
        else if(o.omh) omh.reset(new OrderMinHash(*o.omh));
        else if(o.cfss) cfss.reset(new VSetSketch(*o.cfss));
    }
    OptSketcher(const Dashing2Options &opts): enc_(opts.enc_), enc128_(enc_.to_u128()), ence_(enc_.to_entmin64()), ence128_(enc_.to_entmin128()), k_(opts.k_), w_(opts.w_), use128_(opts.use128()) {
        input_mode(opts.input_mode());
        assert(opts.hashtype() == it_);
        if(opts.sspace_ == SPACE_EDIT_DISTANCE)
            omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
        else if(opts.sketch_compressed_set) {
            switch(int(opts.fd_level_ * 2.)) {
                case 1: cfss.reset(new VSetSketch(NibbleSetS(opts.count_threshold_, opts.sketchsize_, opts.compressed_b_, opts.compressed_a_))); break;
                case 2: cfss.reset(new VSetSketch(ByteSetS(opts.count_threshold_, opts.sketchsize_, opts.compressed_b_, opts.compressed_a_))); break;
                case 4: cfss.reset(new VSetSketch(ShortSetS(opts.count_threshold_, opts.sketchsize_, opts.compressed_b_, opts.compressed_a_))); break;
                case 8: cfss.reset(new VSetSketch(UintSetS(opts.count_threshold_, opts.sketchsize_, opts.compressed_b_, opts.compressed_a_))); break;
                default: []() __attribute__((cold,noinline)) {THROW_EXCEPTION(std::invalid_argument("Unexpected fd_level."));}();
            }
        }
    }
    template<typename Func>
    void for_each(const Func &func, const char *s, size_t n) {
        if(use128_) {
            if(unsigned(k_) <= enc_.nremperres128()) {
                if(entmin)
                    ence128_.for_each(func, s, n);
                else
                    enc128_.for_each(func, s, n);
            } else {
                rh128_.for_each(func, s, n);
            }
        } else {
            if(unsigned(k_) <= enc_.nremperres64()) {
                if(entmin)
                    ence_.for_each(func, s, n);
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

FastxSketchingResult &fastx2sketch_byseq(FastxSketchingResult &ret, Dashing2Options &opts, const std::string &path, kseq_t *kseqs, std::string outpath, bool parallel, size_t seqs_per_batch) {
    gzFile ifp;
    kseq_t *myseq = kseqs ? &kseqs[OMP_ELSE(omp_get_thread_num(), 0)]: (kseq_t *)std::calloc(sizeof(kseq_t), 1);
    size_t batch_index = 0;
    OptSketcher sketcher(opts);
    sketcher.rh_ = opts.rh_;
    sketcher.rh128_ = opts.rh128_;
    if(verbosity >= INFO) {
        std::fprintf(stderr, "Beginning sketch-by-seq with sketch space %s\n", to_string(opts.sspace_).data());
    }
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm()) {
            sketcher.opss.reset(new OPSetSketch(opts.sketchsize_));
            if(opts.count_threshold_ > 0) sketcher.opss->set_mincount(opts.count_threshold_);
        } else if(!opts.sketch_compressed_set)
            sketcher.fss.reset(new FullSetSketch(opts.count_threshold_, opts.sketchsize_, opts.save_kmers_, opts.save_kmercounts_));
        else {
            assert(sketcher.cfss);
        }
    } else if(opts.sspace_ == SPACE_MULTISET) {
        DBG_ONLY(std::fprintf(stderr, "Using BMH\n"););
        sketcher.bmh.reset(new BagMinHash(opts.sketchsize_, opts.save_kmers_, opts.save_kmercounts_));
    } else if(opts.sspace_ == SPACE_PSET) {
        sketcher.pmh.reset(new ProbMinHash(opts.sketchsize_));
        DBG_ONLY(std::fprintf(stderr, "Setting sketcher.pmh: %p\n", (void *)sketcher.pmh.get()););
    } else if(opts.sspace_ == SPACE_EDIT_DISTANCE) {
        sketcher.omh.reset(new OrderMinHash(opts.sketchsize_, opts.k_));
        DBG_ONLY(std::fprintf(stderr, "Setting sketcher.omh: %p\n", (void *)sketcher.omh.get()););
    } else THROW_EXCEPTION(std::runtime_error("Should have been set space, multiset, probset, or edit distance"));
    if(opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_PSET) {
        sketcher.ctr.reset(new Counter(opts.cssize()));
    }
    size_t total_nseqs = 0;
    for_each_substr([&total_nseqs](const std::string &path) {
        gzFile fp = gzopen(path.data(), "r");
        if(!fp) THROW_EXCEPTION(std::runtime_error("Failed to open gzfile"s + path + "to count sequences."));
        kseq_t *ks = kseq_init(fp);
        while(kseq_read(ks) >= 0) ++total_nseqs;
        kseq_destroy(ks);
        gzclose(fp);
    }, path);
    if(verbosity >= INFO) {
        std::fprintf(stderr, "%zu seqs in file\n", total_nseqs);
    }
    if(outpath.size() && outpath != "-" && outpath != "/dev/stdout") {
        if(!bns::isfile(outpath)) {
            DBG_ONLY(std::fprintf(stderr, "Creating outpath '%s'\n", outpath.data()););
            std::FILE *fp = std::fopen(outpath.data(), "wb");
            if(!fp) THROW_EXCEPTION(std::runtime_error("Failed to open path "s + outpath + " for writing"));
            std::fclose(fp);
        }
        size_t offset = 16 + sizeof(double) * total_nseqs;
        if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
            offset = 20 + sizeof(double) * total_nseqs;
        }
        if(int rc = ::truncate(outpath.data(), offset); rc) {
            std::fprintf(stderr, "is outpath %s a file? %d\n", outpath.data(), bns::isfile(outpath));
            THROW_EXCEPTION(std::runtime_error("Failed to resize signature file for fastx2sketch_byseq. rc: "s + std::to_string(rc)));
        }
        ret.signatures_.assign(outpath);
    }

    const size_t num_sigs = (total_nseqs * opts.sketchsize_) >> opts.sigshift();
    if(verbosity >= INFO) {
        std::fprintf(stderr, "Reserving signatures space: %zu entries, for %zu bytes each\n", num_sigs, (sizeof(double) * opts.sketchsize_) >> opts.sigshift());
    }
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE) {
        ret.signatures_.reserve(num_sigs);
    }

    ret.cardinalities_.resize(total_nseqs);
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE && opts.save_kmers_) {
        ret.kmers_.resize(opts.sketchsize_ * total_nseqs);
    }

    if((opts.kmer_result_ != FULL_MMER_SEQUENCE) && opts.save_kmercounts_) {
        ret.kmercounts_.resize(opts.sketchsize_ * total_nseqs);
    }
    if(verbosity >= INFO) {
        std::fprintf(stderr, "Reserved kmers, kmer counts, etc.\n");
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
    DBG_ONLY(std::fprintf(stderr, "save ids: %d, save counts %d\n", opts.save_kmers_, opts.save_kmercounts_););
    size_t lastindex = 0;
    for_each_substr([&](const auto &x) {
        DBG_ONLY(std::fprintf(stderr, "Processing substr %s\n", x.data()););
        if((ifp = gzopen(x.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to read from ") + x));
        gzbuffer(ifp, 1u << 17);
        kseq_assign(myseq, ifp);
        for(int c;(c = kseq_read(myseq)) >= 0;) {
            //DBG_ONLY(std::fprintf(stderr, "Sequence %s of length %zu\n", myseq->name.s, myseq->seq.l););
            ret.sequences_.emplace_back(myseq->seq.s, myseq->seq.l);
            const int off = (myseq->name.s[0] == '>');
            ret.names_.emplace_back(myseq->name.s + off, myseq->name.l - off);
            if(++batch_index == seqs_per_batch) {
                DBG_ONLY(std::fprintf(stderr, "batch index = %zu\n", batch_index););
                resize_fill(opts, ret, std::min(size_t(seqs_per_batch), size_t(total_nseqs - ret.names_.size())), sketching_data, lastindex, nt);
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
    std::fprintf(stderr, "ret signatures size %zu\n", ret.signatures_.size());
    std::fprintf(stderr, "ret signatures capacity %zu\n", ret.signatures_.capacity());
    std::fprintf(stderr, "ret names size %zu\n", ret.names_.size());
#endif
    return ret;
}
void resize_fill(Dashing2Options &opts, FastxSketchingResult &ret, size_t newsz, std::vector<OptSketcher> &sketchvec, size_t &lastindex, size_t nt) {
    if(verbosity >= DEBUG) {
        std::fprintf(stderr, "Calling resize_fill with newsz = %zu\n", newsz);
    }
    const size_t oldsz = ret.names_.size();
    newsz = oldsz + newsz;
    const int sigshift = opts.sigshift();
    if(opts.kmer_result_ != FULL_MMER_SEQUENCE) {
        DBG_ONLY(std::fprintf(stderr, "old sig size %zu, cap %zu, new %zu\n", ret.signatures_.size(), ret.signatures_.capacity(), newsz * opts.sketchsize_ >> sigshift););
        assert(oldsz * (opts.sketchsize_ >> sigshift) <= ret.signatures_.capacity());
        ret.signatures_.resize(oldsz * (opts.sketchsize_ >> sigshift));
        DBG_ONLY(std::fprintf(stderr, "ret signature size: %zu\n", ret.signatures_.size()););
    }
    DBG_ONLY(std::fprintf(stderr, "mmer matrix size %zu. save kmers %d\n", ret.kmers_.size(), opts.save_kmers_););
    DBG_ONLY(std::fprintf(stderr, "Parsing %s\n", sketchvec.front().enable_protein() ? "Protein": "DNA"););
    std::unique_ptr<std::vector<uint64_t>[]> seqmins;
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
        seqmins.reset(new std::vector<uint64_t>[(oldsz - lastindex)]);
    }
    OMP_PRAGMA("omp parallel for num_threads(nt) schedule(dynamic)")
    for(size_t i = lastindex; i < oldsz; ++i) {
        const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        //DBG_ONLY(std::fprintf(stderr, "%zu/%zu -- parsing sequence from tid = %d\n", i, oldsz, tid););
        auto &sketchers(sketchvec[tid]);
        sketchers.reset();
        const auto seqp = ret.sequences_[i].data();
        const auto seql = ret.sequences_[i].size();
        if(sketchers.omh) { // OrderMinHash
            std::vector<uint64_t> res = sketchers.omh->hash(seqp, seql);
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
            ret.cardinalities_[i] = seql;
        } else if(seqmins) { // Sequence of minimizers
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
            }, seqp, seql);
            // This handles the case where the sequence is shorter than the window size
            // and the entire sequence yields no minimizer.
            // We instead take the minimum-hashed value from the b-tree
            if(opts.w_ > opts.k_ && myseq.empty() && seql) {
                //std::fprintf(stderr, "Adding in single item for small seq]\n");
                if(sketchers.rh128_.n_in_queue()) {
                    u128_t v = sketchers.rh128_.max_in_queue().el_;
                    myseq.push_back(0); myseq.push_back(0);
                    std::memcpy(&myseq[myseq.size() - 2], &v, sizeof(v));
                } else if(sketchers.enc128_.n_in_queue()) {
                    u128_t v = sketchers.enc128_.max_in_queue().el_;
                    myseq.push_back(0); myseq.push_back(0);
                    std::memcpy(&myseq[myseq.size() - 2], &v, sizeof(v));
                } else if(sketchers.ence128_.n_in_queue()) {
                    u128_t v = sketchers.ence128_.max_in_queue().el_;
                    myseq.push_back(0); myseq.push_back(0);
                    std::memcpy(&myseq[myseq.size() - 2], &v, sizeof(v));
                } else if(sketchers.rh_.n_in_queue()) {
                    myseq.push_back(sketchers.rh_.max_in_queue().el_);
                } else if(sketchers.enc_.n_in_queue()) {
                    myseq.push_back(sketchers.enc_.max_in_queue().el_);
                } else if(sketchers.ence_.n_in_queue()) {
                    myseq.push_back(sketchers.ence_.max_in_queue().el_);
                }
            }
            ret.cardinalities_[i] = myseq.size();
        } else {
            DBG_ONLY(std::fprintf(stderr, "Sketching all k-mers in sequence %s/%zu\n", ret.sequences_[i].data(), i););
            const bool isop = sketchers.opss.get(), isctr = sketchers.ctr.get(), isfs = sketchers.fss.get(), iscfss = sketchers.cfss.get();
            auto fsfunc = [&](auto x) __attribute__((always_inline)) {
                x = maskfn(x);
                if(opts.fs_->in_set(x)) return;
                if(isop)    sketchers.opss->update(x);
                else if(isctr) sketchers.ctr->add(x);
                else if(isfs) sketchers.fss->update(x);
                else if(iscfss) std::visit([&x](auto &sketch) __attribute__((always_inline)) {sketch.update(x);}, *sketchers.cfss);
            };
            auto nofsfunc = [&](auto x) __attribute__((always_inline)) {
                x = maskfn(x);
                if(isop) sketchers.opss->update(x);
                else if(isctr) sketchers.ctr->add(x);
                else if(isfs) sketchers.fss->update(x);
                else if(iscfss) std::visit([&x](auto &sketch) __attribute__((always_inline)) {sketch.update(x);}, *sketchers.cfss);
            };
            if(opts.fs_) {
                sketchers.for_each(fsfunc, seqp, seql);
            } else {
                sketchers.for_each(nofsfunc, seqp, seql);
            }
            DBG_ONLY(std::fprintf(stderr, "Sketched all k-mers in sequence\n"););
            RegT *ptr = nullptr;
            const uint64_t *kmer_ptr = nullptr;
            std::vector<double> kmercounts;
            if(opts.sspace_ == SPACE_SET) {
                if(sketchers.ctr) {
                    DBG_ONLY(std::fprintf(stderr, "Sketching setspace with count thresholding.\n"););
                    if(sketchers.opss) {
                        sketchers.ctr->finalize(*sketchers.opss, opts.count_threshold_);
                        ptr = sketchers.opss->data();
                        ret.cardinalities_[i] = sketchers.opss->getcard();
                    } else if(sketchers.fss) {
                        sketchers.ctr->finalize(*sketchers.fss, opts.count_threshold_);
                        ptr = sketchers.fss->data();
                        ret.cardinalities_[i] = sketchers.fss->getcard();
                    } else if(sketchers.cfss) {
                        std::visit([&](auto &sketch) {sketchers.ctr->finalize(sketch, opts.count_threshold_); ptr = (RegT *)sketch.data();ret.cardinalities_[i] = sketch.getcard();}, *sketchers.cfss);
                    }
                    DBG_ONLY(std::fprintf(stderr, "Sketched setspace with count thresholding.\n"););
                } else {
                    ptr = sketchers.opss ? sketchers.opss->data(): sketchers.fss ? sketchers.fss->data(): (RegT *)getdata(*sketchers.cfss);
                    ret.cardinalities_[i] = sketchers.opss ? sketchers.opss->getcard(): sketchers.fss ? sketchers.fss->getcard(): sketchers.cfss ? getcard(*sketchers.cfss): std::numeric_limits<double>::quiet_NaN();
                    if(std::isnan(ret.cardinalities_[i])) {
                        if(verbosity >= INFO) {
                            std::fprintf(stderr, "Warning: sketch cardinality was a NAN. This is unexpected. Settings to 0.\n");
                        }
                        ret.cardinalities_[i] = 0.;
                    }
                    if(ret.cardinalities_[i] < 10 * opts.sketchsize_) {
                        DBG_ONLY(std::fprintf(stderr, "Cardinality exact counting fall-back\n"););
                        flat_hash_set<uint64_t> ids;
                        ids.reserve(opts.sketchsize_);
                        if(opts.fs_) {
                            sketchers.for_each([&](auto x) {
                                x = maskfn(x);
                                if(opts.fs_->in_set(x)) return;
                                ids.insert(x);
                            }, seqp, seql);
                        } else {
                            sketchers.for_each([&](auto x) {ids.insert(maskfn(x));}, seqp, seql);
                        }
                        ret.cardinalities_[i] = ids.size();
                        DBG_ONLY(std::fprintf(stderr, "Cardinality exact counting fall-back complete\n"););
                    }
                    kmer_ptr = sketchers.opss ? sketchers.opss->ids().data(): sketchers.fss ? sketchers.fss->ids().data(): (const uint64_t *)nullptr;
                    if(sketchers.opss && sketchers.opss->idcounts().size()) {
                        auto &idc = sketchers.opss->idcounts();
                        kmercounts.resize(opts.sketchsize_);
                        std::copy(idc.begin(), idc.end(), kmercounts.begin());
                    } else if(sketchers.fss && sketchers.fss->idcounts().size()) {
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
            if(!sketchers.cfss || opts.fd_level_ != 0.5) {
                const size_t mss = opts.sketchsize_ * i;
                DBG_ONLY(std::fprintf(stderr, "About to copy out. For sketch idx %zu, Offset into dict: %zu. copying %zu bytes.\n", mss, mss >> sigshift, ((opts.sketchsize_ * sizeof(RegT)) >> sigshift)););
                std::memcpy(&ret.signatures_[mss >> sigshift], ptr, ((opts.sketchsize_ * sizeof(RegT)) >> sigshift));
                // std::copy(ptr, ptr + (opts.sketchsize_ >> sigshift), &ret.signatures_[i * (opts.sketchsize_ >> sigshift)]);
                DBG_ONLY(std::fprintf(stderr, "Copied out. For sketch idx %zu, Offset into dict: %zu. copying %zu bytes.\n", mss, mss >> sigshift, ((opts.sketchsize_ * sizeof(RegT)) >> sigshift)););
            } else {
                const uint8_t *const srcptr = std::get<NibbleSetS>(*sketchers.cfss).data();
                uint8_t *destptr = (uint8_t *)&ret.signatures_[i * opts.sketchsize_ >> sigshift];
                for(size_t i = 0; i < opts.sketchsize_; i += 2) {
                    *destptr++ = (srcptr[i] << 4) | srcptr[i + 1];
                }
            }
            if(kmer_ptr && ret.kmers_.size()) {
                std::copy(kmer_ptr, kmer_ptr + opts.sketchsize_, &ret.kmers_[i * opts.sketchsize_]);
            } else if(verbosity >= DEBUG) {
                std::fprintf(stderr, "k-mers not saved\n");
            }
            if(kmercounts.size() && ret.kmercounts_.size()) {
                std::copy(kmercounts.begin(), kmercounts.end(), &ret.kmercounts_[i * opts.sketchsize_]);
            }
        }
    }
    if(seqmins) {
        const size_t seqminsz = oldsz - lastindex;
        using OT = std::conditional_t<(sizeof(RegT) == 8), uint64_t, std::conditional_t<(sizeof(RegT) == 4), uint32_t, u128_t>>;
        static_assert(sizeof(RegT) == sizeof(OT), "Size of hash registers must match that being sketched");
        size_t new_sigsz = std::accumulate(seqmins.get(), seqmins.get() + seqminsz, size_t(0), [](size_t cs, const auto &v) {return cs + v.size();});
        if constexpr(sizeof(RegT) == 16) new_sigsz >>= 1;
        else if constexpr(sizeof(RegT) == 4) new_sigsz <<= 1;
        std::vector<size_t> offsets{{ret.signatures_.size()}};
        ret.signatures_.resize(new_sigsz + ret.signatures_.size());
        for(size_t i = 0; i < seqminsz; ++i) {
            auto noff = seqmins[i].size();
            if constexpr(sizeof(RegT) == 16) {
                noff >>= 1;
            }
            offsets.push_back(offsets.back() + noff);
        }
        size_t oldnpf = ret.nperfile_.size();
        ret.nperfile_.resize(oldnpf + seqminsz);
        OMP_PFOR
        for(size_t i = 0; i < seqminsz; ++i) {
            auto &x(seqmins[i]);
            if(opts.use128_) {
                auto ptr = (u128_t *)x.data();
                for(size_t K = 0; K < (x.size() >> 1); ++K) {
                    u128_t item;
                    std::memcpy(&item, ptr + K, sizeof(item));
                    item = invmaskfn(item);
                    std::memcpy(ptr + K, &item, sizeof(item));
                }
            } else {
                std::transform(x.begin(), x.end(), x.begin(), [](auto item) {return invmaskfn(item);});
            }
            const OT *ptr = (const OT *)x.data();
            size_t xsz = x.size();
            if constexpr(sizeof(RegT) == 16) {
                xsz >>= 1;
            }
            ret.nperfile_[oldnpf + i] = xsz;
            std::memcpy(&ret.signatures_[offsets[i]], ptr, sizeof(RegT) * xsz);
        }
    }
    lastindex = oldsz;
}
} // dashing2
