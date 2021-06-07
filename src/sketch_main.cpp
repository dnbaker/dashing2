#include "sketch_core.h"
#include "options.h"
#include "cmp_main.h"



#define SKETCH_OPTS \
static option_struct sketch_long_options[] = {\
    LO_FLAG("canon", 'C', canon, true)\
    LO_FLAG("cache", 'W', cache, true)\
    LO_FLAG("multiset", OPTARG_DUMMY, s, SPACE_MULTISET)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
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
    LO_ARG("save-kmers", 's')\
    LO_ARG("ffile", 'F')\
    LO_ARG("qfile", 'Q')\
    LO_ARG("count-threshold", 'm')\
    LO_ARG("outfile", 'o')\
    SHARED_OPTS\
};



namespace dashing2 {
void sketch_usage() {
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
                         "-2/--128bit/long-kmers: Use 128-bit k-mer hashes instead of 64-bit\n"
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
                         "--prob: Sketch m-mers into ProbMinHash. Treats weighted sets as discrete probability distributions.\n"
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
                         //"-G/--seq: Full m-mer sequence. This faster than building the hash set, and can be used to build a minimizer index afterwards\n"
                         //"          On the other hand, it can require higher memory for large sequence collections\n"
                         "\nMetadata Options\n"
                         "If sketching, you can also choose to save k-mers (the IDs corresponding to the k-mer selected), or\n"
                         " and optionally save the counts for these k-mers\n"
                         "This could be used to build inverted indexes (using samples to estimate containment), or for frequency estimation\n"
                         "-s/--save-kmers: Save m-mers. This puts the m-mers saved into .kmer files to correspond with the minhash samples.\n"
                         "-N/--save-kmercounts: Save m-mer counts for sketches. This puts the m-mer counts saved into .kmercounts.f64 files to correspond with the m-mers.\n"
                         "-o/--outfile: sketches are stacked into a single file and written to [arg]\n"
                         SHARED_DOC_LINES
    );
}


int sketch_main(int argc, char **argv) {
    int c;
    int k = 16, w = -1, nt = -1;
    SketchSpace s = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = false;
    bool exact_kmer_dist = false;
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
    size_t nbytes_for_fastdists = sizeof(RegT);
    bool parse_by_seq = false;
    Measure measure = SIMILARITY;
    // By default, use full hash values, but allow people to enable smaller
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    SKETCH_OPTS
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:Ns2BPWh?ZJGH", sketch_long_options, &option_index)) >= 0;) {
        switch(c) {
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
        case OPTARG_ISZ: measure = INTERSECTION; break;
        case OPTARG_OUTPREF: {
            outprefix = optarg; break;
        }
        SHARED_FIELDS
        case '?': case 'h': sketch_usage(); return 1;
    }}
    if(nt < 0) {
        if(char *s = std::getenv("OMP_NUM_THREADS")) nt = std::atoi(s);
        else nt = 1;
    }
    std::fprintf(stderr, "rest: %s\n", to_string(res).data());
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
        .save_kmercounts(save_kmercounts)
        .parse_by_seq(parse_by_seq);
    std::fprintf(stderr, "opts save kmers: %d\n", opts.save_kmers_);
    opts.count_threshold_ = count_threshold;
    if(opts.sspace_ == SPACE_PSET && opts.kmer_result_ == ONE_PERM) opts.kmer_result_ = FULL_SETSKETCH;
    if(paths.empty()) {
        std::fprintf(stderr, "No paths provided. See usage.\n");
        sketch_usage();
        return 1;
    }
    SketchingResult result = sketch_core(opts, paths, outfile);
    result.nqueries(nq); // TODO: use nqueries to perform asymmetric comparisons
    if(cmpout.size()) {
        Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist);
        distopts.measure_ = measure;
        cmp_core(distopts, result);
    }
    return 0;
}

} // namespace dashing2
