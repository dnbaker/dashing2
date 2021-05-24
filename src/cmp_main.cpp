#include "cmp_main.h"
#include "sketch_core.h"
#include "options.h"
#include "refine.h"

namespace dashing2 {

#define CMP_OPTS(name) \
static option_struct name[] = {\
    LO_FLAG("canon", 'C', canon, true)\
    LO_FLAG("cache", 'W', cache, true)\
    LO_FLAG("multiset", OPTARG_DUMMY, s, SPACE_MULTISET)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, s, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', s, SPACE_PSET)\
    LO_FLAG("edit-distance", 'E', s, SPACE_EDIT_DISTANCE)\
    LO_FLAG("set", 'H', res, FULL_MMER_SET)\
    LO_FLAG("full-setsketch", 'Z', res, FULL_SETSKETCH)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
    LO_FLAG("save-kmers", 's', save_kmers, true)\
    LO_FLAG("enable-protein", OPTARG1, rht, bns::RollingHashingType::PROTEIN)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
    LO_FLAG("presketched", OPTARG_PRESKETCHED, presketched, true)\
    LO_ARG("regbytes", OPTARG_REGBYTES)\
    LO_ARG("outprefix", OPTARG_OUTPREF)\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("window-size", 'w')\
    LO_ARG("countsketch-size", 'c')\
    LO_ARG("threads", 'p')\
    LO_ARG("sketchsize", 'S')\
    LO_ARG("save-kmercounts", 'N')\
    LO_ARG("ffile", 'F')\
    LO_ARG("qfile", 'Q')\
    LO_ARG("count-threshold", 'm')\
    LO_ARG("outfile", 'o')\
    SHARED_OPTS\
}

void cmp_usage() {
    std::fprintf(stderr, "dashing2 cmp usage is not written.\n");
}
void load_results(Dashing2DistOptions &opts, SketchingResult &result, const std::vector<std::string> &paths) {
    auto &pf = paths.front();
    struct stat st;
    ::stat(pf.data(), &st);
    if(paths.size() == 1) {
        std::string namesf = pf + ".names.txt";
        if(bns::isfile(namesf)) {
            std::string l;
            for(std::ifstream ifs(namesf);std::getline(ifs, l);) {
                result.names_.emplace_back(l);
            }
        } else throw std::runtime_error("cmp expects a packed sketch matrix or multiple paths");
        assert(result.names_.size());
        //size_t nregs = st.st_size / result.names_.size() / sizeof(RegT);
        std::FILE *fp = std::fopen(pf.data(), "w");
        if(!fp) throw std::runtime_error(std::string("Failed to open ") + pf);
        result.signatures_.resize(st.st_size / sizeof(RegT));
        if(std::fread(result.signatures_.data(), st.st_size, 1, fp) != 1u) throw std::runtime_error(std::string("Failed to write block of size ") + std::to_string(st.st_size));
        std::fclose(fp);
    } else { // Else, we have to load sketches from each file
        result.nperfile_.resize(paths.size());
        auto &fsizes = result.nperfile_;
        std::vector<size_t> csizes(fsizes.size());
        for(size_t i = 0; i < paths.size(); ++i) {
            struct stat st;
            ::stat(paths[i].data(), &st);
            assert(st.st_size % sizeof(RegT) == 0);
            fsizes[i] = st.st_size / sizeof(RegT);
            csizes[i + 1] = csizes[i] + st.st_size;
        }
        const size_t totalsize = std::accumulate(fsizes.begin(), fsizes.end(), size_t(0));
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
                if((ifp = std::fopen(p.data(), "rb")) == nullptr) throw std::runtime_error(std::string("Could not open path ") + p);
                if(std::fread(&result.kmers_[csizes[i]], sizeof(uint64_t), fsizes[i], ifp) != fsizes[i]) {
                    std::fprintf(stderr, "Failed to read at path %s\n", p.data());
                    std::fclose(ifp);
                    std::exit(1);
                }
                std::fclose(ifp);
            }
            if(std::string p(path + ".kmercounts.f64"); bns::isfile(p)) {
                if((ifp = std::fopen(p.data(), "rb")) == nullptr) throw 2;
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

#if 0
enum KmerSketchResultType {
    ONE_PERM = 0,       // Faster (3-4x) than Full, comparable accuracy for both cardinality and set similarities
    FULL_SETSKETCH = 1, // Not stochastically-averaged; potentially better LSH properties
    FULL_MMER_SET = 2,
    FULL_MMER_SEQUENCE = 3,
    FULL_MMER_COUNTDICT = 4
};

#endif

int cmp_main(int argc, char **argv) {
    int c;
    int k = 16, w = 50, nt = -1;
    SketchSpace s = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = false, presketched = false;
    double count_threshold = 0., similarity_threshold = -1.;
    size_t cssize = 0, sketchsize = 1024;
    std::string ffile, outfile, qfile;
    int option_index = 0;
    bns::RollingHashingType rht = bns::DNA;
    DataType dt = DataType::FASTX;
    std::string outprefix;
    OutputKind ok = SYMMETRIC_ALL_PAIRS;
    std::string cmpout; // Only used if distances are also requested
    int topk_threshold = -1;
    int truncate_mode = 0;
    double nbytes_for_fastdists = sizeof(RegT);
    bool parse_by_seq = false;
    // By default, use full hash values, but allow people to enable smaller
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    CMP_OPTS(cmp_long_options);
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:Ns2BPWh?ZJGH", cmp_long_options, &option_index)) >= 0;) {switch(c) {
        case 'k': k = std::atoi(optarg); break;
        case 'w': w = std::atoi(optarg); break;
        case 'W': cache = true; break;
        case 'B': s = SPACE_MULTISET; res = FULL_SETSKETCH; break;
        case 'P': s = SPACE_PSET; res = FULL_SETSKETCH; break;
        case 'Z': res = FULL_SETSKETCH; break;
        case 'o': outfile = optarg; break;
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;
        case 'C': canon = true; break;
        case 'p': nt = std::atoi(optarg); break;
        case 'S': sketchsize = std::atoi(optarg); break;
        case 'N': save_kmers = save_kmercounts = true; break;
        case 's': save_kmers = true; break;
        case 'H': res = FULL_MMER_SET; break;
        case 'J': res = FULL_MMER_COUNTDICT; break;
        case 'G': res = FULL_MMER_SEQUENCE; break;
        case '2': use128 = true; break;
        case 'm': count_threshold = std::atof(optarg); break;
        case 'F': ffile = optarg; break;
        case 'Q': qfile = optarg; break;
        case OPTARG_REGBYTES: nbytes_for_fastdists = std::atof(optarg); break;
        case OPTARG_OUTPREF: {
            outprefix = optarg; break;
        }
        SHARED_FIELDS
        case '?': case 'h': cmp_usage(); return 1;
    }}
    if(nt < 0) {
        if(char *s = std::getenv("OMP_NUM_THREADS")) nt = std::atoi(s);
        else nt = 1;
    }
    std::vector<std::string> paths(argv + optind, argv + argc);
    std::unique_ptr<std::vector<std::string>> qup;
    if(ffile.size()) {
        std::ifstream ifs(ffile);
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
    Dashing2Options opts(k, w, rht, s, dt);
    if(nbytes_for_fastdists == 0.5 && (opts.sketchsize_ % 2)) {
        std::fprintf(stderr, "Increasing sketch size to a multiple of 2 to avoid cutting in half\n");
        ++opts.sketchsize_;
    }
    opts.nthreads(nt)
        .kmer_result(res)
        .cache_sketches(cache)
        .cssize(cssize)
        .use128(use128)
        .sketchsize(sketchsize)
        .save_kmers(save_kmers)
        .outprefix(outprefix)
        .save_kmercounts(save_kmercounts)
        .parse_by_seq(parse_by_seq);
    opts.count_threshold_ = count_threshold;
    Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout);
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
    }
    cmp_core(distopts, result);
    return 0;
}


}
