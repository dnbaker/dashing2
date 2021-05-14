#include <getopt.h>
#include "fastxsketch.h"
#include "bwsketch.h"
#include "bedsketch.h"
#include "lfsketch.h"


#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

using option_struct = struct option;

enum OptArg{
    OPTARG1 = 1000,
    OPTARG_BED,
    OPTARG_BIGWIG,
    OPTARG_LEAFCUTTER,
    OPTARG_DUMMY
};

#define SKETCH_OPTS \
static option_struct sketch_long_options[] = {\
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
    LO_FLAG("save-kmers", 's', save_kmers, true)\
    LO_FLAG("enable-protein", OPTARG1, rht, bns::RollingHashingType::PROTEIN)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
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
};


using namespace dashing2;
size_t nbytes_from_line(const std::string &line) {
    size_t ret = 0;
    for_each_substr([&ret](const std::string &s) {ret += bns::filesize(s.data());}, line);
    return ret;
}


void usage() {
    std::fprintf(stderr, "dashing2 sketch <opts> [fastas... (optional)]\n"
                         "We use only m-mers; if w <= k, however, this reduces to k-mers\n"
                         "Flags:\n\n"
                         "Runtime options:\n"
                         "-p/--threads: Set number of threads [1]\n"
                         "Encoding options\n"
                         "Dashing2 can sketch 4 kinds of files:"
                         "Fastq/Fasta, which has specific encoding options (default)\n"
                         "--bed to sketch BED files for interval sets\n"
                         "--bigwig to sketch BigWig files for coverage vectors\n"
                         "and --leafcutter to sketch LeafCutter splicing output\n"
                         "\n\nFastx Options:\n"
                         "-k/--kmer-length: set k\n"
                         "-w/--window-size: set window size for winnowing; by default, all m-mers are used.\n"
                         "-2: Use 128-bit k-mer hashes instead of 64-bit\n"
                         "-m/--threshold: Set a count threshold for inclusion. Default: 0.\n"
                         "--enable-protein: Switch from DNA-sequence encoding to protein encoding. This treats all characters as valid\n"
                         "\nPathsOptions\n\n"
                         "By default, dashing2 reads positional arguments and sketches them. You may want to use flags instructing it\n"
                         "to read from paths in <file>. Additionally, you can put multiple files separated by spaces into a single line "
                         "to place them all into a single sketch instead.\n"
                         "-F/--ffile: read paths from file in addition to positional arguments\n"
                         "-Q/--qfile: read query paths from file; this is used for asymmetric queries (e.g., containment)\n"
                         "This accelerates weighted sketching at the cost of some approximation.\n"
                         "\nSketch options\n"
                         "These decide how m-mers are accumulated.\n"
                         "Default behavior is set sketching (tossing multiplicities). If --multiset or --prob is set or a minimum count is provided,"
                         "\nk-mers will be counted before sketching.\n"
                         "-S/--sketchsize: Set sketchsize (1024)\n"
                         "In sketching space you can use ProbMinHash, BagMinHash, or SetSketch, which is set MinHash\n"
                         "-P/--prob: Sketch m-mers into ProbMinHash. Treats weighted sets as discrete probability distributions.\n"
                         "-B/--multiset: Sketch m-mers into BagMinHash. Treats weighted sets as multisets.\n"
                         "-Z/--full-setsketch: Full setsketch (not stochastically-averaged)\n"
                         "This should perform similarly to default setsketch behavior, but has better behaviors with large sketches and small sets\n"
                         "It typically comes at 2-4x runtime cost, depending on sketch size\n"
                         "-c/--countsketch-size: Use Count-Sketch counting instead of exact counting, using [arg] as the size.\n    "
                         "This allows you to avoid unbounded dictionary size at the cost of some approximation of weighted sets\n"
                         "This only affects methods which perform counting\n"
                         "You can also emit full m-mer sets, a count dictionary (key-count map)\n"
                         "-H/--set: Full m-mer set. This generates a sorted hash set for m-mers in the data. If the parser is windowed (-w is set), this may be rather small.\n"
                         "-J/--countdict: Full m-mer countdict. This generates a sorted hash set for m-mers in the data, and additionally saves the associated counts for these m-mers.\n"
                         "-G/--seq: Full m-mer sequence. This faster than building the hash set, and can be used to build a minimizer index afterwards\n"
                         "          On the other hand, it can require higher memory for large sequence collections\n"
                         "\nMetadata Options\n"
                         "If sketching, you can also choose to save k-mers (the IDs corresponding to the k-mer selected), or\n"
                         " and optionally save the counts for these k-mers\n"
                         "This could be used to build inverted indexes (using samples to estimate containment), or for frequency estimation\n"
                         "-s/--save-kmers: Save m-mers. This puts the m-mers saved into .kmer files to correspond with the minhash samples.\n"
                         "-N/--save-kmercounts: Save m-mer counts for sketches. This puts the m-mer counts saved into .kmercounts.f32 files to correspond with the m-mers.\n"
                         "-o/--outfile: sketches are stacked into a single file and written to [arg]\n"
    );
}

int sketch_main(int argc, char **argv) {
    int c;
    int k = 16, w = 50, nt = -1;
    SketchSpace s = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = false;
    double threshold = 0.;
    size_t cssize = 0, sketchsize = 1024;
    std::string ffile, outfile, qfile;
    int option_index = 0;
    bns::RollingHashingType rht = bns::DNA;
    DataType dt = DataType::FASTX;
    std::string outprefix;
    SKETCH_OPTS
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:Ns2BPWh?ZJGH", sketch_long_options, &option_index)) >= 0;) {switch(c) {
        case 'k': k = std::atoi(optarg); break;
        case 'w': w = std::atoi(optarg); break;
        //case 'W': cache = true; break;
        //case 'B': s = SPACE_MULTISET; res = FULL_SETSKETCH; break;
        //case 'P': s = SPACE_PSET; res = FULL_SETSKETCH; break;
        //case 'Z': res = FULL_SETSKETCH; break;
        case 'o': outfile = optarg; break;
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;
        //case 'C': canon = true; break;
        case 'p': nt = std::atoi(optarg); break;
        case 'S': sketchsize = std::atoi(optarg); break;
        case 'N': save_kmers = save_kmercounts = true; break;
        //case 's': save_kmers = true; break;
        //case 'H': res = FULL_MMER_SET; break;
        //case 'J': res = FULL_MMER_COUNTDICT; break;
        //case 'G': res = FULL_MMER_SEQUENCE; break;
        //case '2': use128 = true; break;
        case 'm': threshold = std::atof(optarg); break;
        case 'F': ffile = optarg; break;
        case 'Q': qfile = optarg; break;
        case OPTARG_OUTPREF: {
            outprefix = optarg; break;
        }
        case '?': case 'h': usage(); return 1;
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
    std::fprintf(stderr, "Sketching %zu arguments (lhs) and %zu (rhs)\n", nref, nq);
    Dashing2Options opts(k, w, rht, s, dt);
    opts.nthreads(nt)
        .kmer_result(res)
        .cache_sketches(cache)
        .cssize(cssize)
        .use128(use128)
        .sketchsize(sketchsize)
        .save_kmers(save_kmers)
        .outprefix(outprefix)
        .save_kmercounts(save_kmercounts);
    opts.count_threshold_ = threshold;
    if(paths.empty()) {
        std::fprintf(stderr, "No paths provided. See usage.\n");
        usage();
        return 1;
    }
    const size_t npaths = paths.size();
    std::vector<std::pair<size_t, uint32_t>> filesizes(npaths);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        filesizes[i] = {nbytes_from_line(paths[i]), uint32_t(i)};
    }
    std::sort(filesizes.begin(), filesizes.end(), std::greater<>());
    SketchingResult result;
    if(opts.dtype_ == DataType::FASTX) {
        result = fastx2sketch(opts, paths);
    } else if(opts.dtype_ == DataType::LEAFCUTTER) {
        auto res = lf2sketch(paths, opts);
        result.signatures_ = std::move(res.registers());
        result.names_ = std::move(res.sample_names());
        result.nperfile_.resize(res.nsamples_per_file().size());
        std::copy(res.nsamples_per_file().begin(), res.nsamples_per_file().end(), result.nperfile_.begin());
    } else if(opts.dtype_ == DataType::BED || opts.dtype_ == DataType::BIGWIG) {
        result.signatures_.resize(npaths * opts.sketchsize_);
        OMP_PFOR_DYN
        for(size_t i = 0; i < npaths; ++i) {
            auto myind = filesizes.size() ? filesizes[i].second: uint32_t(i);
            std::vector<RegT> sigs;
            if(opts.dtype_ == DataType::BED) {
                sigs = bed2sketch(paths[myind], opts);
                std::copy(sigs.begin(), sigs.end(), &result.signatures_[myind * opts.sketchsize_]);
            } else if(opts.dtype_ == DataType::BIGWIG) {
                if(opts.by_chrom_) {
                    std::fprintf(stderr, "Warning: by_chrom is ignored for bigwig sketching. Currently, all sets are grouped together. To group by chromosome, split the BW file by chromosome.");
                    opts.by_chrom_ = false;
                }
                auto res = bw2sketch(paths[myind], opts);
                sigs = std::move(*res.global_.get());
            }
            std::copy(sigs.begin(), sigs.end(), &result.signatures_[myind * opts.sketchsize_]);
        }
    }
    if(paths.size() == 1 && outfile.empty()) {
        const std::string suf = 
                opts.sspace_ == SPACE_SET ? (opts.one_perm_ ? ".opss": ".ss"):
                opts.sspace_ == SPACE_MULTISET ? ".bmh":
                opts.sspace_ == SPACE_PSET ? ".pmh" :
                opts.sspace_ == SPACE_EDIT_DISTANCE ? ".omh": ".unknown_sketch";
        outfile = paths.front() + suf;
        if(opts.trim_folder_paths_) {
            outfile = trim_folder(path);
            if(opts.outprefix_)
                outfile = opts.outprefix_ + '/' + outfile;
        }
    }
    if(outfile.size()) {
        if(result.signatures_.empty()) throw std::runtime_error("Can't write stacked sketches if signatures were not generated");
        std::fprintf(stderr, "Writing stacked sketches to %s\n", outfile.data());
        std::FILE *ofp = std::fopen(outfile.data(), "wb");
        if(!ofp) throw std::runtime_error(std::string("Failed to open file at ") + outfile);
        std::fwrite(result.signatures_.data(), sizeof(RegT), result.signatures_.size(), ofp);
        std::fclose(ofp);
        if(result.names_.size()) {
            if((ofp = std::fopen((outfile + ".names.txt").data(), "wb")) == nullptr)
                throw std::runtime_error(std::string("Failed to open outfile at ") + outfile + ".names.txt");
            for(const auto &n: result.names_) {
                if(std::fwrite(n.data(), 1, n.size(), ofp) != n.size()) throw std::runtime_error("Failed to write names to file");
                std::fputc('\n', ofp);
            }
            std::fclose(ofp);
        }
    }
#ifndef NDEBUG
    if(result.names_.size()) {
        for(size_t i = 0; i < result.names_.size(); ++i) {
            std::fprintf(stderr, "%zu/%zu: %s\n", i + 1, result.names_.size(), result.names_[i].data());
        }
    }
#endif
    return 0;
}
