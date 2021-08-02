#include "cmp_main.h"
#include "sketch_core.h"
#include "options.h"
#include "refine.h"
#include <filesystem>

namespace dashing2 {

#define CMP_OPTS(name) \
static option_struct name[] = {\
    LO_FLAG("presketched", OPTARG_PRESKETCHED, presketched, true)\
    SHARED_OPTS\
}



void cmp_usage() {
    std::fprintf(stderr, "dashing2 cmp usage is not written.\n");
}
void load_results(Dashing2DistOptions &opts, SketchingResult &result, const std::vector<std::string> &paths) {
    std::fprintf(stderr, "Loading results using Dashing2Options: %s\n", opts.to_string().data());
    auto &pf = paths.front();
    struct stat st;
    ::stat(pf.data(), &st);
    if(paths.size() == 1) {
        std::string namesf = pf + ".names.txt";
        if(bns::isfile(namesf)) {
            std::string l;
            for(std::ifstream ifs(namesf);std::getline(ifs, l);) {
                std::cerr << "Line: " << l << '\n';
                if(l.empty() || l.front() == '#') continue;
                const typename std::string::size_type it = l.find_first_of('\t');
                result.names_.emplace_back(l.substr(0, it));
                if(it != std::string::npos) {
                    char *s = &l[it + 1];
                    result.cardinalities_.emplace_back(std::strtod(s, &s));
                    if(*s) {
                        std::fprintf(stderr, "Saving kmer count file %s\n", s + 1);
                        result.kmercountfiles_.push_back(s + 1);
                    }
                }
            }
        } else {
            THROW_EXCEPTION(std::runtime_error("cmp expects a packed sketch matrix or multiple paths"););
        }
        assert(result.names_.size());
        //size_t nregs = st.st_size / result.names_.size() / sizeof(RegT);
        std::FILE *fp = std::fopen(pf.data(), "w");
        if(!fp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open ") + pf));
        result.signatures_.resize(st.st_size / sizeof(RegT));
        if(std::fread(result.signatures_.data(), st.st_size, 1, fp) != 1u) THROW_EXCEPTION(std::runtime_error(std::string("Failed to write block of size ") + std::to_string(st.st_size)));
        std::fclose(fp);
    } else { // Else, we have to load sketches from each file
        result.nperfile_.resize(paths.size());
        auto &fsizes = result.nperfile_;
        std::vector<size_t> csizes(fsizes.size() + 1);
        for(size_t i = 0; i < paths.size(); ++i) {
            struct stat st;
            ::stat(paths[i].data(), &st);
            assert(st.st_size % sizeof(RegT) == 0);
            const auto nelem = st.st_size / sizeof(RegT);
            fsizes[i] = nelem;
            csizes[i + 1] = csizes[i] + nelem;
        }
        const size_t totalsize = csizes.back();
        result.signatures_.resize(totalsize); // Account for the size of the sketch registers
        if(bns::isfile(pf + ".kmerhashes.u64")) {
            std::fprintf(stderr, "Loading k-mer hashes, too\n");
            result.kmers_.resize(result.signatures_.size());
        }
        if(bns::isfile(pf + ".kmercounts.f64")) {
            std::fprintf(stderr, "Loading k-mer counts, too\n");
            result.kmercounts_.resize(result.signatures_.size());
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < paths.size(); ++i) {
            auto &path = paths[i];
            std::FILE *ifp = std::fopen(path.data(), "rb");
            if(std::fread(&result.signatures_[csizes[i]], sizeof(RegT), fsizes[i], ifp) != fsizes[i]) {
                std::fprintf(stderr, "Failed to read at path %s\n", path.data());
                std::fclose(ifp);
                std::exit(1);
            }
            std::fclose(ifp);
            if(std::string p(path + ".kmerhashes.u64"); bns::isfile(p)) {
                if((ifp = std::fopen(p.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Could not open path ") + p));
                if(std::fread(&result.kmers_[csizes[i]], sizeof(uint64_t), fsizes[i], ifp) != fsizes[i]) {
                    std::fprintf(stderr, "Failed to read at path %s\n", p.data());
                    std::fclose(ifp);
                    std::exit(1);
                }
                std::fclose(ifp);
            }
            if(std::string p(path + ".kmercounts.f64"); bns::isfile(p)) {
                if((ifp = std::fopen(p.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open kmer count file"));
                if(std::fread(&result.kmercounts_[csizes[i]], sizeof(double), fsizes[i], ifp) != fsizes[i]) {
                    std::fprintf(stderr, "Failed to read k-mer counts at path %s\n", p.data());
                    std::fclose(ifp);
                    std::exit(1);
                }
                std::fclose(ifp);
            }
        }
    }
    std::fprintf(stderr, "loading results, but this isn't written!\n");
}

int cmp_main(int argc, char **argv) {
    int c;
    int k = 16, w = 0, nt = -1;
    SketchSpace sketch_space = SPACE_SET;
    KmerSketchResultType res = FULL_SETSKETCH;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = true, presketched = false;
    bool exact_kmer_dist = false;
    bool refine_exact = false; // This uses sketching for K-NN graph generation, then uses exact distances for NN refinement
    unsigned int count_threshold = 0.;
    double similarity_threshold = -1.;
    size_t cssize = 0, sketchsize = 1024;
    std::string ffile, outfile, qfile;
    int option_index = 0;
    bns::RollingHashingType rht = bns::DNA;
    DataType dt = DataType::FASTX;
    std::string outprefix;
    OutputKind ok = SYMMETRIC_ALL_PAIRS;
    std::string cmpout; // Only used if distances are also requested
    bool normalize_bed = false;
    int topk_threshold = -1;
    int truncate_mode = 0;
    double nbytes_for_fastdists = sizeof(RegT);
    double downsample_frac = 1.;
    bool parse_by_seq = false;
    bool hpcompress = false;
    std::string fsarg;
    Measure measure = SIMILARITY;
    uint64_t seedseed = 0;
    size_t batch_size = 16;
    std::string spacing;
    // By default, use full hash values, but allow people to enable smaller
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    CMP_OPTS(cmp_long_options);
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:Ns2BPWh?ZJGH", cmp_long_options, &option_index)) >= 0;) {switch(c) {
        SHARED_FIELDS
        case OPTARG_HELP: case '?': case 'h': cmp_usage(); return 1;
    }}
    std::vector<std::string> paths(argv + optind, argv + argc);
    std::unique_ptr<std::vector<std::string>> qup;
    std::string cmd(std::filesystem::absolute(std::filesystem::path(argv[-1])));
    for(char **s = argv; *s; cmd += std::string(" ") + *s++);
    std::fprintf(stderr, "[Dashing2] Invocation: %s ", cmd.data());
    if(nt < 0) {
        char *s = std::getenv("OMP_NUM_THREADS");
        if(s) nt = std::max(std::atoi(s), 1);
    }
    OMP_ONLY(omp_set_num_threads(nt));
    if(ffile.size()) {
        std::ifstream ifs(ffile);
        static constexpr size_t bufsize = 1<<18;
        std::unique_ptr<char []> buf(new char[bufsize]);
        ifs.rdbuf()->pubsetbuf(buf.get(), bufsize);
        for(std::string l;std::getline(ifs, l);) {
            paths.push_back(l);
        }
    }
    size_t nref = paths.size();
    if(qfile.size()) {
        std::ifstream ifs(qfile);
        for(std::string l;std::getline(ifs, l);)
            paths.push_back(l);
    }
    size_t nq = paths.size() - nref;
    Dashing2Options opts(k, w, rht, sketch_space, dt, nt, use128, spacing, canon, res);
    opts
        .cache_sketches(cache)
        .cssize(cssize)
        .sketchsize(sketchsize)
        .save_kmers(save_kmers)
        .outprefix(outprefix)
        .save_kmercounts(save_kmercounts)
        .save_kmers(save_kmers)
        .parse_by_seq(parse_by_seq)
        .cmd(cmd).count_threshold(count_threshold)
        .homopolymer_compress_minimizers(hpcompress)
        .seedseed(seedseed);
    if(hpcompress) {
        if(!opts.homopolymer_compress_minimizers_) THROW_EXCEPTION(std::runtime_error("Failed to hpcompress minimizers"));
    }
    opts.filterset(fsarg);
    if(nbytes_for_fastdists == 0.5)
        opts.sketchsize_ += opts.sketchsize_ & 1; // Ensure that sketch size is a multiple of 2 if using nibbles
    opts.bed_parse_normalize_intervals_ = normalize_bed;
    opts.downsample(downsample_frac);
    Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist, refine_exact);
    distopts.measure_ = measure;
    distopts.cmp_batch_size_ = std::max(batch_size, size_t(distopts.nthreads()));
    SketchingResult result;
    if(presketched) {
        std::set<std::string> suffixset;
        for(const auto &p: paths) {
            suffixset.insert(p.substr(p.find_last_of('.'), std::string::npos));
        }
        std::string suf;
        if(suffixset.size() != 1) {
            std::fprintf(stderr, "Multiple suffixes in set (%s). Picking the first one.\n", suffixset.begin()->data());
            suf = paths.front().substr(paths.front().find_last_of('.'), std::string::npos);
        } else suf = *suffixset.begin();
        if(suf == ".bmh") {
            distopts.sspace_ = SPACE_MULTISET;
            distopts.kmer_result(FULL_SETSKETCH);
        } else if(suf == ".pmh") {
            distopts.sspace_ = SPACE_PSET;
            distopts.kmer_result(FULL_SETSKETCH);
        } else if(suf == ".ss" || suf == ".opss") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(suf == ".opss" ? ONE_PERM: FULL_SETSKETCH);
        } else if(suf == ".kmerset64") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(FULL_MMER_SET);
            distopts.use128(false);
        } else if(suf == ".kmerset128") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(FULL_MMER_SET);
            distopts.use128(true);
        } else if(suf == "mmerseq64") {
            distopts.sspace_ = SPACE_SET;
            auto &path = paths.front();
            std::string countg = path.substr(0, path.find_last_of('.')) + "kmercounts.f64";
            if(bns::isfile(countg)) {
                distopts.kmer_result(FULL_MMER_COUNTDICT);
            } else distopts.kmer_result(FULL_MMER_SEQUENCE);
            distopts.use128(false);
        } else if(suf == "mmerseq128") {
            distopts.sspace_ = SPACE_SET;
            auto &path = paths.front();
            std::string countg = path.substr(0, path.find_last_of('.')) + "kmercounts.f64";
            if(bns::isfile(countg)) {
                distopts.kmer_result(FULL_MMER_COUNTDICT);
            } else distopts.kmer_result(FULL_MMER_SEQUENCE);
            distopts.use128(true);
        }
        load_results(distopts, result, paths);
    } else {
        result = sketch_core(distopts, paths, outfile);
        result.nqueries(nq);
    }
    cmp_core(distopts, result);
    return 0;
}


}
