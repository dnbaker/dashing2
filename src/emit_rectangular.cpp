#include "cmp_main.h"

namespace dashing2 {

void emit_rectangular(Dashing2DistOptions &opts, const SketchingResult &result) {
    const size_t ns = result.names_.size();
    std::FILE *ofp = opts.outfile_path_.empty() || opts.outfile_path_ == "-" ? stdout: std::fopen(opts.outfile_path_.data(), "w");
    if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open path at ") + opts.outfile_path_));
    std::setvbuf(ofp, nullptr, _IOFBF, 1<<17);
    const bool asym = opts.output_kind_ == ASYMMETRIC_ALL_PAIRS;
    std::deque<std::pair<std::unique_ptr<float[]>, size_t>> datq;
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
        std::fprintf(ofp, "%zu\n", ns);
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
            } else {
                auto &f = datq.front();
                const size_t nwritten = opts.output_kind_ == PANEL ? nq: asym ? ns: ns - f.second - 1;
                if(opts.output_format_ == HUMAN_READABLE) {
                    auto datp = f.first.get();
                    auto &seqn = result.names_[f.second];
                    std::string fn(seqn);
                    if(fn.size() < 9) fn.append(9 - fn.size(), ' ');
                    std::fwrite(fn.data(), 1, fn.size(), ofp);
                    for(size_t i = 0; i < nwritten;std::fprintf(ofp, "\t%0.9g", datp[i++]));
                    std::fputc('\n', ofp);
                } else if(opts.output_format_ == MACHINE_READABLE) {
                    if(std::fwrite(datq.front().first.get(), sizeof(float), nwritten, ofp) != nwritten)
                        THROW_EXCEPTION(std::runtime_error(std::string("Failed to write row ") + std::to_string(datq.front().second) + " to disk"));
                } // else throw std::runtime_error("This should never happen");
                std::lock_guard<std::mutex> guard(datq_lock);
                datq.pop_front();
            }
        }
    });
    if(opts.output_kind_ == PANEL) {
        for(size_t i = 0; i < nf; ++i) {
            std::unique_ptr<float[]> dat(new float[nq]);
            OMP_PFOR_DYN
            for(size_t j = 0; j < nq; ++j) {
                dat[j] = compare(opts, result, i, j + nf);
            }
            std::lock_guard<std::mutex> guard(datq_lock);
            datq.emplace_back(std::pair<std::unique_ptr<float[]>, size_t>{std::move(dat), i});
        }
    } else {
        for(size_t i = 0; i < ns; ++i) {
            // TODO: batch queries together for cache efficiency (distmat::parallel_fill for an example)
            size_t nelem = asym ? ns: ns - i - 1;
            std::unique_ptr<float[]> dat(new float[nelem]);
            const auto datp = dat.get() - (asym ? size_t(0): i + 1);
            OMP_PFOR_DYN
            for(size_t start = asym ? 0: i + 1;start < ns; ++start) {
                datp[start] = compare(opts, result, i, start);
            }
            std::lock_guard<std::mutex> guard(datq_lock);
            datq.emplace_back(std::pair<std::unique_ptr<float[]>, size_t>{std::move(dat), i});
        }
    }
    loopint = 1;
    if(sub.joinable()) sub.join();
    if(ofp != stdout) std::fclose(ofp);
}


} // namespace dashing2
