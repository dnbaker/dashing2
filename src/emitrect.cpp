#include "cmp_main.h"
#include "fmt/format.h"
#include "fmt/os.h"
#include <optional>

namespace dashing2 {
using namespace std::literals::string_literals;

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
    const auto &nwritten() const {return std::get<3>(*this);}
};

template<size_t L>
constexpr std::array<char, 2 * L + 1> make_tablut() {
    std::array<char, 2 * L + 1> ret{0};
    for(size_t i = 0; i < L; ++i) {
        size_t id = i * 2;
        ret[id] = '\t';
        ret[id + 1] = '-';
    }
    ret[2 * L] = '\0';
    return ret;
}

static constexpr const std::array<char, 513> tabarr = make_tablut<256>();
static constexpr std::string_view tabstr(tabarr.data(), 512);

static INLINE int print_tabs(size_t n, std::FILE *ofp) {
    static constexpr const char *ts = tabarr.data();
    while(n > 256) {
        if(std::fwrite(ts, 512, 1, ofp) != 1u) return -1;
        n -= 256;
    }
    const auto nw = 2u * n;
    return std::fwrite(ts, 1, nw, ofp) == nw ? 0: -2;
}

static INLINE void print_tabs(size_t n, std::back_insert_iterator<fmt::memory_buffer> &biof) {
    for(;n > 256; format_to(biof, tabstr), n -= 256);
    const auto substr = tabstr.substr(0, n << 1);
    format_to(biof, substr);
}

static INLINE void print_tabs(size_t n, fmt::ostream &os) {
    for(;n > 256; os.print("{}", tabstr), n -= 256);
    const auto substr = tabstr.substr(0, n << 1);
#ifndef NDEBUG
    for(size_t i = 0; i < substr.size() >> 1; ++i) {
        assert(substr[i * 2] == '\t');
        assert(substr[i * 2 + 1] == '-');
    }
#endif
    os.print("{}", substr);
}
#ifndef NDEBUG
#define EMPTY -137.
#endif
using cfloatp = const float *;
#if USE_FMT_FPRINTF
static constexpr const char *FMT8 = "\t{:0.7g}\t{:0.7g}\t{:0.7g}\t{:0.7g}\t{:0.7g}\t{:0.7g}\t{:0.7g}\t{:0.7g}";
static constexpr const char *FMT1 = "\t{:0.7g}";
#else
static constexpr const char *FMT8 = "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}";
static constexpr const char *FMT1 = "\t{}";
#endif

void batched_write(cfloatp &src, fmt::ostream &of, const size_t jend) {
    // Batched formatting provides a significant speed advantage
    const float *dat8end = src + (jend / 8) * 8;
    const float *datend = src + jend;
    for(;src < dat8end; src += 8)
        of.print(FMT8,
                 *src, src[1], src[2], src[3], src[4], src[5], src[6], src[7]);

    for(;src < datend; ++src)
        of.print(FMT1, *src);
    of.print("\n");
}

void batched_write(const float * &src, std::back_insert_iterator<fmt::memory_buffer> &biof, const size_t jend) {
    // Batched formatting provides a significant speed advantage
    const float *dat8end = src + (jend / 8) * 8;
    const float *datend = src + jend;
    for(;src < dat8end; src += 8)
        fmt::format_to(biof, FMT8,
            *src, src[1], src[2], src[3], src[4], src[5], src[6], src[7]);

    for(;src < datend; ++src) {
        fmt::format_to(biof, FMT1, *src);
    }
    fmt::format_to(biof, "\n");
}

void emit_rectangular(const Dashing2DistOptions &opts, const SketchingResult &result) {
    const size_t ns = result.names_.size();
    const std::string outp = opts.outfile_path_.empty() || opts.outfile_path_.front() == '-'
            ? "/dev/stdout"s: opts.outfile_path_;
    // Only make fmt::ostream if emitting in human-readable form
    std::optional<fmt::ostream> ofopt(opts.output_format_ == HUMAN_READABLE
                ? std::optional<fmt::ostream>(fmt::output_file(outp, fmt::buffer_size=131072))
                : std::optional<fmt::ostream>(std::nullopt));
    std::FILE *ofp = 0;
    if(opts.output_format_ == MACHINE_READABLE) {
        if(opts.outfile_path_.empty() || opts.outfile_path_.front() == '-') {
            ofp = stdout;
            buffer_to_blksize(ofp);
        } else {
            if((ofp = bfopen(opts.outfile_path_.data(), "wb")) == 0)
                THROW_EXCEPTION(std::runtime_error("Failed to open path "s + opts.outfile_path_ + " for writing"));
        }
    }
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::deque<QTup> datq;
    volatile int loopint = 0;
    std::mutex datq_lock;
#ifndef NDEBUG
    if(opts.output_format_ == MACHINE_READABLE) {
        std::fprintf(stderr, "Emitting machine-readable: %s\n", to_string(opts.output_format_).data());
    } else {
        std::fprintf(stderr, "Emitting human-readable: %s\n", to_string(opts.output_format_).data());
    }
#endif
    // Emit Header
    if(opts.output_format_ == HUMAN_READABLE) {
        auto &of = ofopt.value();
        if(opts.output_kind_ != PHYLIP) {
            const char *labelstr = asym ? "Asymmetric pairwise": opts.output_kind_ == PANEL ? "Panel (Query/Refernce)": "Symmetric pairwise";
            of.print("#Dashing2 {} Output\n", labelstr);
            of.print("#Dashing2Options: {}\n", opts.to_string());
            of.print("#Sources");
            for(size_t i = 0; i < result.names_.size(); of.print("\t{}", result.names_[i++]));
            of.print("\n");
        } else {
            of.print("{}\n", ns);
        }
    }
    /*
     *  This is a worker thread which processes and emits
     *  data in the datq which computation is done in parallel by other threads.
     *  If the queue is empty, it sleeps.
     */
    const size_t nq = result.nqueries(), nf = ns - nq;
    std::thread sub = std::thread([&](){
        while(loopint == 0) {
            if(datq.empty()) {
                std::this_thread::sleep_for(std::chrono::duration<size_t, std::nano>(500));
                continue;
            }
            auto &f = datq.front();
            auto fs = f.start(), fe = f.stop();
#ifndef NDEBUG
            if(loopint) {
                std::fprintf(stderr, "Writing data from queue of size %zu after the loopint is terminated. This means all computation is done and we are just formatting and emitting data.\n", datq.size());
            }
#endif
            if(opts.output_format_ == HUMAN_READABLE) {
                auto &of = ofopt.value();
                const float *datp = f.data();
                for(size_t i = fs; i < fe; ++i) {
                    auto &seqn = result.names_[i];
                    std::string fn(seqn);
                    if(fn.size() < 9) fn.append(9 - fn.size(), ' ');
                    of.print("{}", fn);
                    const size_t jend = opts.output_kind_ == PANEL ? nq: asym ? ns: ns - i - 1;
                    if(opts.output_kind_ == SYMMETRIC_ALL_PAIRS) print_tabs(i + 1, of);
                    batched_write(datp, of, jend);
                }
            } else {
                assert(opts.output_format_ == MACHINE_READABLE);
                const size_t nwritten = datq.front().nwritten();
                if(std::fwrite(datq.front().data(), sizeof(float), nwritten, ofp) != nwritten)
                    THROW_EXCEPTION(std::runtime_error(std::string("Failed to write rows ") + std::to_string(datq.front().start()) + "-" + std::to_string(datq.front().stop()) + " to disk"));
            }
            std::lock_guard<std::mutex> guard(datq_lock);
            datq.pop_front();
        }
    });
    const size_t batch_size = std::max(std::min(unsigned(opts.cmp_batch_size_), opts.nthreads()), 1u);
    // We have two access patterns --
    // Unbatched (batch_size <= 1), which fills in the matrix one row at a time
    // and
    // Batched (batch_size > 1), which is more cache efficient by grouping comparisons
    // so that computations using the same data can share it
    if(opts.output_kind_ == PANEL) {
        if(batch_size <= 1) {
            for(size_t i = 0; i < nf; ++i) {
                std::unique_ptr<float[]> dat(new float[nq]);
#ifndef NDEBUG
                std::fill_n(dat.get(), nq, EMPTY);
#endif
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
#ifndef NDEBUG
                std::fill_n(dat.get(), nwritten, EMPTY);
#endif
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
                const size_t diff = erow - firstrow;
                const size_t nwritten = ns * diff;
                std::unique_ptr<float[]> dat(new float[nwritten]);
#ifndef NDEBUG
                std::fill_n(dat.get(), nwritten, EMPTY);
#endif
                OMP_PFOR_DYN
                for(size_t fs = firstrow; fs < erow; ++fs) {
                    auto datp = &dat[(fs - firstrow) * ns];
                    for(size_t j = 0; j < ns; ++j)
                        datp[j] = compare(opts, result, fs, j);
                }
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.emplace_back(QTup{std::move(dat), firstrow, erow, nwritten});
            }
        } else { // all-pairs symmetric! (upper-triangular)
            if(batch_size <= 1) {
                for(size_t i = 0; i < ns; ++i) {
                    size_t nelem = asym ? ns: ns - i - 1;
                    std::unique_ptr<float[]> dat(new float[nelem]);
#ifndef NDEBUG
                    std::fill_n(dat.get(), nelem, EMPTY);
#endif
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
                for(size_t bi = 0; bi < nbatches; ++bi) {
                    const size_t firstrow = bi * batch_size;
                    const size_t erow = std::min((bi + 1) * batch_size, ns);
                    std::vector<size_t> offsets{0};
                    DBG_ONLY(size_t sum = 0;)
                    for(size_t fs = firstrow; fs < erow; ++fs) {
                        offsets.push_back(ns - fs - 1 + offsets.back());
                        DBG_ONLY(sum += (ns - fs - 1);)
                    }
                    const size_t nwritten = offsets.back();
                    assert(nwritten == sum);
                    auto dat = std::make_unique<float[]>(nwritten);
                    DBG_ONLY(std::fill_n(dat.get(), nwritten, EMPTY);)
                    std::atomic<size_t> totalused(0);
                    OMP_PFOR_DYN
                    for(size_t fs = firstrow; fs < erow; ++fs) {
                        auto myoff = fs - firstrow;
#ifndef NDEBUG
                        size_t shouldoff = 0;
                        for(size_t ofs = firstrow; ofs < fs; ++ofs) shouldoff += ns - ofs - 1;
                        assert(shouldoff == offsets.at(myoff) || !std::fprintf(stderr, "Expected %zu for offsets, found %zu\n", shouldoff, offsets.at(myoff)));
#endif
                        auto datp = &dat[offsets[myoff]] - fs - 1;
                        for(size_t j = fs + 1; j < ns; ++j) {
                            DBG_ONLY(++totalused;)
                            datp[j] = compare(opts, result, fs, j);
                        }
                    }
                    assert(totalused.load() == nwritten);
                    assert(std::all_of(dat.get(), dat.get() + nwritten, [](auto x) {return !std::isinf(x);}));
                    std::lock_guard<std::mutex> guard(datq_lock);
                    datq.emplace_back(QTup{std::move(dat), firstrow, erow, nwritten});
                }
            }
        }
    }
    loopint = 1;
    if(sub.joinable()) sub.join();
    if(const size_t dqs = datq.size(); dqs) {
        if(opts.output_format_ == HUMAN_READABLE) {
            auto &of = ofopt.value();
            of.flush();
            std::vector<fmt::memory_buffer> bufs(dqs);
            // Format in parallel, then write in sequence
            OMP_PFOR_DYN
            for(size_t di = 0; di < dqs; ++di) {
                auto &f = datq[di];
                auto biof = std::back_inserter(bufs[di]);
                const auto fs = f.start(), fe = f.stop();
                const float *datp = f.data();
                for(size_t i = fs; i < fe; ++i) {
                    std::string fn(result.names_[i]);
                    if(fn.size() < 9) fn.append(9 - fn.size(), ' ');
                    format_to(biof, fn);
                    const size_t jend = opts.output_kind_ == PANEL ? nq: asym ? ns: ns - i - 1;
                    if(opts.output_kind_ == SYMMETRIC_ALL_PAIRS) print_tabs(i + 1, biof);
                    batched_write(datp, biof, jend);
                }
            }
            datq.clear();
            of.close(); // Flush and close, and then open a posix file
            auto off = fmt::file(outp, fmt::file::WRONLY | fmt::file::APPEND);
            for(const auto &buf: bufs) {
                const size_t nblocks = (buf.size() + 131071) / 131072;
                // Chunk the writing according to blocksize
                for(size_t i = 0; i < nblocks; ++i) {
                    const size_t n_in_block = (i == nblocks - 1) ? buf.size() & size_t(131071): size_t(131072);
                    off.write(buf.data() + i * 131072, n_in_block);
                }
            }
        } else {
#define MANUAL_BUFFER 1
#if MANUAL_BUFFER
            const int fd = ::fileno(ofp);
            std::fflush(ofp);
#endif
            while(!datq.empty()) {
                const auto &pair = datq.front();
                const size_t nwritten = pair.nwritten();
#if !MANUAL_BUFFER
                if(std::fwrite(pair.data(), sizeof(float), nwritten, ofp) != nwritten)
                   THROW_EXCEPTION(std::runtime_error(std::string("Failed to write rows ") + std::to_string(pair.start()) + "-" + std::to_string(pair.stop()) + " to disk"));
#else
                const size_t nbytes = sizeof(float) * nwritten;
                static constexpr size_t nfloats_per_block = 131072 / sizeof(float);
                const ssize_t nblocks = (nbytes + 131071) >> 17;
                for(ssize_t i = 0; i < nblocks - 1; ++i) {
                    const ssize_t wrc = ::write(fd, pair.data() + i * nfloats_per_block, 131072ul);
                    if(unlikely(wrc < ssize_t(131072)))
                        throw std::runtime_error("Failed to POSIX write. Wrote "s + std::to_string(wrc) + " instead of 131072");
                }
                const ssize_t lon = (nwritten & 32767ul) * sizeof(float);
                const ssize_t wrc = ::write(fd, pair.data() + (nblocks - 1) * nfloats_per_block, lon);
                if(wrc != lon) {
                   THROW_EXCEPTION(std::runtime_error(std::string("Failed to write rows ") + std::to_string(pair.start()) + "-" + std::to_string(pair.stop()) + " to disk"));
                }
#endif
                datq.pop_front();
            }
        }
    }
    assert(datq.empty());
    if(ofp && ofp != stdout) std::fclose(ofp);
}


} // namespace dashing2
