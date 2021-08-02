#include "cmp_main.h"

namespace dashing2 {

struct QTup: public std::tuple<std::unique_ptr<float[]>, size_t, size_t, size_t> {
    using super = std::tuple<std::unique_ptr<float[]>, size_t, size_t, size_t>;
    template<typename...Args>
    QTup(Args &&...args): super(std::forward<Args>(args)...) {}
    auto data() {return std::get<0>(*this).get();}
    auto data() const {return std::get<0>(*this).get();}
    auto &ptr() {return std::get<0>(*this);}
    auto &ptr() const {return std::get<0>(*this);}
    auto &start() {return std::get<1>(*this);}
    const auto start() const {return std::get<1>(*this);}
    auto &stop() {return std::get<2>(*this);}
    const auto stop() const {return std::get<2>(*this);}
    auto &nwritten() {return std::get<3>(*this);}
    const auto nwritten() const {return std::get<3>(*this);}
};

void emit_rectangular(Dashing2DistOptions &opts, const SketchingResult &result) {
    const size_t ns = result.names_.size();
    std::FILE *ofp = opts.outfile_path_.empty() || opts.outfile_path_ == "-" ? stdout: std::fopen(opts.outfile_path_.data(), "w");
    if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open path at ") + opts.outfile_path_));
    std::setvbuf(ofp, nullptr, _IOFBF, 1<<17);
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::deque<QTup> datq;
    volatile int loopint = 0;
    std::mutex datq_lock;
    if(opts.output_format_ == MACHINE_READABLE) {
        std::fprintf(stderr, "Emitting machine-readable: %s\n", to_string(opts.output_format_).data());
    } else {
        std::fprintf(stderr, "Emitting human-readable: %s\n", to_string(opts.output_format_).data());
    }
    // Emit Header
    if(opts.output_format_ == HUMAN_READABLE) {
        std::fprintf(ofp, "#Dashing2 %s Output\n", asym ? "Asymmetric pairwise": "PHYLIP pairwise");
        std::fprintf(ofp, "#Dashing2Options: %s\n", opts.to_string().data());
        std::fprintf(ofp, "#Sources");
        for(size_t i = 0; i < result.names_.size(); ++i) {
            std::fprintf(ofp, "\t%s", result.names_[i].data());
        }
        std::fputc('\n', ofp);
        if(opts.output_kind_ == SYMMETRIC_ALL_PAIRS) std::fprintf(ofp, "%zu\n", ns);
    }
    /*
     *  This is a worker thread which processes and emits
     *  data in the datq which computation is done in parallel by other threads.
     *  If the queue is empty, it sleeps.
     */
    const size_t nq = result.nqueries(), nf = ns - nq;
    std::thread sub = std::thread([&](){
        while(loopint == 0 || datq.size()) {
            if(datq.empty()) {
                std::this_thread::sleep_for(std::chrono::duration<size_t, std::nano>(2000));
                continue;
            }
            auto &f = datq.front();
            auto fs = f.start(), fe = f.stop();
            if(opts.output_format_ == HUMAN_READABLE) {
                auto datp = f.data();
                for(size_t i = fs; i < fe; ++i) {
                    auto &seqn = result.names_[i];
                    std::string fn(seqn);
                    if(fn.size() < 9) fn.append(9 - fn.size(), ' ');
                    std::fwrite(fn.data(), 1, fn.size(), ofp);
                    const size_t jend = opts.output_kind_ == PANEL ? nq: asym ? ns: ns - i - 1;
                    for(size_t j = 0; j < jend; ++j)
                        std::fprintf(ofp, "\t%0.9g", *datp++);
                    std::fputc('\n', ofp);
                }
            } else if(opts.output_format_ == MACHINE_READABLE) {
                const size_t nwritten = datq.front().nwritten();
                if(std::fwrite(datq.front().data(), sizeof(float), nwritten, ofp) != nwritten)
                    THROW_EXCEPTION(std::runtime_error(std::string("Failed to write rows ") + std::to_string(datq.front().start()) + "-" + std::to_string(datq.front().stop()) + " to disk"));
            }
            std::lock_guard<std::mutex> guard(datq_lock);
            datq.pop_front();
        }
    });
    const size_t batch_size = std::max(opts.cmp_batch_size_, size_t(1));
    if(opts.output_kind_ == PANEL) {
        if(batch_size <= 1) {
            for(size_t i = 0; i < nf; ++i) {
                std::unique_ptr<float[]> dat(new float[nq]);
                OMP_PFOR_DYN
                for(size_t j = 0; j < nq; ++j) {
                    dat[j] = compare(opts, result, i, j + nf);
                }
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.emplace_back(QTup{std::move(dat), i, i + 1, nq});
            }
        } else {
            const size_t nbatches = (nf + batch_size - 1) / batch_size;
            for(size_t bi = 0; bi < nbatches; ++bi) {
                const size_t firstrow = bi * batch_size;
                const size_t erow = std::min((bi + 1) * batch_size, nf);
                const size_t nrow = erow - firstrow;
                const size_t nwritten = nq * nrow;
                std::unique_ptr<float[]> dat(new float[nwritten]);
                OMP_PFOR_DYN
                for(size_t i = 0; i < nrow; ++i) {
                    for(size_t j = 0; j < nq; ++j)
                        dat[i * nq + j] = compare(opts, result, i + firstrow, j + nf);
                }
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.emplace_back(QTup{std::move(dat), firstrow, erow, nwritten});
            }
        }
    } else {
        if(asym) {
            const size_t nbatches = (ns + batch_size - 1) / batch_size;
            for(size_t bi = 0; bi < nbatches; ++bi) {
                const size_t firstrow = bi * batch_size;
                const size_t erow = std::min((bi + 1) * batch_size, ns);
                std::fprintf(stderr, "Batch %zu, ranges from %zu to %zu\n", bi, firstrow, erow);
                const size_t diff = erow - firstrow;
                const size_t nwritten = ns * diff;
                std::unique_ptr<float[]> dat(new float[nwritten]);
                OMP_PFOR_DYN
                for(size_t fs = firstrow; fs < erow; ++fs) {
                    auto datp = &dat[(fs - firstrow) * ns];
                    for(size_t j = 0; j < ns; ++j) {
                        datp[j] = compare(opts, result, fs, j);
                    }
                }
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.emplace_back(QTup{std::move(dat), firstrow, erow, nwritten});
            }
        } else { // all-pairs symmetric! (upper-triangular)
            if(batch_size <= 1) {
                for(size_t i = 0; i < ns; ++i) {
                    size_t nelem = asym ? ns: ns - i - 1;
                    std::unique_ptr<float[]> dat(new float[nelem]);
                    const auto datp = dat.get() - (asym ? size_t(0): i + 1);
                    OMP_PFOR_DYN
                    for(size_t start = asym ? 0: i + 1;start < ns; ++start) {
                        datp[start] = compare(opts, result, i, start);
                    }
                    std::lock_guard<std::mutex> guard(datq_lock);
                    datq.emplace_back(QTup{std::move(dat), i, i + 1, nelem});
                }
            } else {
                const size_t nbatches = (ns + batch_size - 1) / batch_size;
                std::fprintf(stderr, "Gothere\n");
                for(size_t bi = 0; bi < nbatches; ++bi) {
                
                    const size_t firstrow = bi * batch_size;
                    const size_t erow = std::min((bi + 1) * batch_size, ns);
                    std::fprintf(stderr, "Batch %zu, ranges from %zu to %zu\n", bi, firstrow, erow);
                    std::vector<size_t> offsets{0};
                    for(size_t fs = firstrow; fs < erow; ++fs) {
                        offsets.push_back(ns - fs - 1 + offsets.back());
                    }
                    const size_t nwritten = std::accumulate(offsets.begin(), offsets.end(), size_t(0));
                    std::fprintf(stderr, "%zu written, for %zu as batch size\n", nwritten, batch_size);
                    std::unique_ptr<float[]> dat(new float[nwritten]);
                    std::fill_n(dat.get(), nwritten, -13.f);
                    OMP_PFOR_DYN
                    for(size_t fs = firstrow; fs < erow; ++fs) {
                        auto myoff = fs - firstrow;
                        std::fprintf(stderr, "row %zu at offset %zu\n", fs, offsets[myoff]);
                        size_t shouldoff = 0;
                        for(size_t ofs = firstrow; ofs < fs; ++ofs) shouldoff += ns - ofs - 1;
                        assert(shouldoff == offsets[myoff] || !std::fprintf(stderr, "Expected %zu for offsets, found %zu\n", shouldoff, offsets[myoff]));
                        auto datp = &dat[offsets[myoff]];
                        for(size_t j = fs; ++j < ns;) {
                            std::fprintf(stderr, "Data at offset %zu is being set fo %zu/%zu\n", j - fs, fs, j);
                            datp[j - fs] = compare(opts, result, fs, j);
                        }
                    }
                    std::lock_guard<std::mutex> guard(datq_lock);
                    datq.emplace_back(QTup{std::move(dat), firstrow, erow, nwritten});
                }
            }
        }
    }
    loopint = 1;
    if(sub.joinable()) sub.join();
    if(ofp != stdout) std::fclose(ofp);
}


} // namespace dashing2
