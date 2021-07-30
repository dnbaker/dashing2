#include "lfsketch.h"

namespace dashing2 {

LFResult LFResult::merge_results(const LFResult *start, size_t n) {
    LFResult ret;
    std::vector<size_t> offsets(n + 1);
    ret.nsamples_per_file().resize(n);
    for(size_t i = 1; i < n + 1; ++i) {
        const size_t sz = start[i - 1].nsamples_per_file().front();
        ret.nsamples_per_file()[i - 1] = sz;
        offsets[i] = offsets[i - 1] + sz;
    }
    const size_t total_samples = offsets.back();
    const size_t ss = start->registers().size() / start->sample_names().size();
    ret.registers().resize(total_samples * ss);
    auto &destnames = ret.sample_names();
    destnames.resize(total_samples);
    // Merge the names together, merge samples together
    OMP_PFOR
    for(size_t i = 0; i < n; ++i) {
        const size_t beg = offsets[i];
        auto &srcfn = start[i].filenames().front();
        size_t pos = srcfn.find("_perind");
        if(pos == std::string::npos) {
            if((pos = srcfn.find(".count")) == std::string::npos)
            pos = srcfn.find("_");
        }
        std::string pref = srcfn.substr(0, pos);
        auto &srcnames = start[i].sample_names();
        std::transform(srcnames.begin(), srcnames.end(), &destnames[beg],
                       [pref](const auto &s) {return s + ':' + pref;});
        std::copy(start[i].registers().begin(), start[i].registers().end(), &ret.registers()[beg]);
    }
    return ret;
}
static constexpr size_t GZ_BUFFER_SIZE = 524288;
LFResult lf2sketch(std::string path, const Dashing2Options &opts) {
    if(opts.sspace_ > SPACE_PSET) THROW_EXCEPTION(std::invalid_argument("Can't do edit distance for Splice junction files"));
    LFResult ret;
    ret.filenames() = {path};
    std::unique_ptr<char[]> gzbuf(new char[GZ_BUFFER_SIZE]);
    gzFile ifp;
    if((ifp = gzopen(path.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open input file ") + path));
    char *line;
    if(!(line = gzgets(ifp, gzbuf.get(), GZ_BUFFER_SIZE))) THROW_EXCEPTION(std::runtime_error("Failed to read line from gzFile... is it empty?"));

    for(char *s = std::strchr(line, ' ') + 1;s;) {
        char *s2 = std::strchr(s, ' ');
        if(!s2) {
            size_t nal = 0;
            for(s2 = s;!std::isspace(*s2) && *s2; ++s2) {
                nal += std::isalnum(*s2);
            }
            if(!nal) break;
        }
        ret.sample_names().emplace_back(s, s2);
        if(*s2 == 0) break;
        s = s2 + 1;
    }
    const size_t nsamples = ret.sample_names().size();
    ret.nsamples_per_file() = {nsamples};
    std::unique_ptr<std::vector<FullSetSketch>> ss;
    std::unique_ptr<std::vector<OPSetSketch>> opss;
    std::unique_ptr<std::vector<BagMinHash>> bmhs;
    std::unique_ptr<std::vector<ProbMinHash>> pmhs;
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm()) {
            opss.reset(new std::vector<OPSetSketch>());
            while(opss->size() < nsamples) opss->emplace_back(opts.sketchsize_);
        } else
            ss.reset(new std::vector<FullSetSketch>(nsamples, FullSetSketch(opts.count_threshold_, opts.sketchsize_)));
    } else if(opts.sspace_ == SPACE_MULTISET) {
        bmhs.reset(new std::vector<BagMinHash>(nsamples, BagMinHash(opts.sketchsize_)));
    } else if(opts.sspace_ == SPACE_PSET) {
        pmhs.reset(new std::vector<ProbMinHash>(nsamples, ProbMinHash(opts.sketchsize_)));
    }
    //auto t = std::chrono::high_resolution_clock::now();
    //size_t ln = 0;
    for(char *line;(line = gzgets(ifp, gzbuf.get(), GZ_BUFFER_SIZE)) != nullptr;) {
        //if(++ln % 1024 == 0) std::fprintf(stderr, "%zu lines read in %gms\n", ln, std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t).count());
        char *lend = line;
        int ncolons = 0;
        while(*lend && ncolons < 3) ncolons += (*lend++ == ':');
        if(opts.trim_chr_ && (*line == 'c' || *line == 'C') && std::memcmp(line + 1, "hr", 2) == 0)
            line += 3;
        std::string splice_site(line, lend - 1);
        const uint64_t splice_hash = std::hash<std::string>{}(splice_site);
        ret.splice_sites().push_back(std::move(splice_site));
        for(size_t sample_id = 0;(lend = std::strchr(lend, ' ')) != nullptr;++sample_id) {
            double num = std::strtoul(lend + 1, &lend, 10);
            if(num == 0) continue;
            unsigned long denom = std::strtoul(lend + 1, &lend, 10);
            if(ss) (*ss)[sample_id].update(splice_hash);
            else if(opss) (*opss)[sample_id].update(splice_hash);
            else {
                if(opts.bed_parse_normalize_intervals_) num /= denom;
                if(bmhs)
                    (*bmhs)[sample_id].update(splice_hash, num);
                else //if(pmhs)
                    (*pmhs)[sample_id].update(splice_hash, num);
            }
        }
    }
    gzclose(ifp);
    ret.registers().resize(nsamples * opts.sketchsize_);
    for(size_t i = 0;i < nsamples; ++i) {
        const RegT *ptr =
            bmhs ? (*bmhs)[i].data():
            pmhs ? (*pmhs)[i].data():
            ss   ? (*ss)[i].data():
            opss   ? (*opss)[i].data(): static_cast<const RegT *>(nullptr);
        if(!ptr) {
            std::fprintf(stderr, "Missing ptr for sample %zu/%s\n", i, ret.sample_names()[i].data());
            continue;
        }
        std::copy(ptr, ptr + opts.sketchsize_, &ret.registers()[opts.sketchsize_ * i]);
    }
    return ret;
}

LFResult lf2sketch(std::vector<std::string> paths, const Dashing2Options &opts) {
    std::vector<LFResult> ret(paths.size());
    OMP_PFOR_DYN
    for(size_t i = 0; i < paths.size(); ++i)
        ret[i] = lf2sketch(paths[i], opts);
    return LFResult::merge_results(ret.data(), ret.size());
}


} // namespace dashing2
