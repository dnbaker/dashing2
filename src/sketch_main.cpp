#include "sketch_core.h"
#include "options.h"
#include "cmp_main.h"
#include <filesystem>



#define SKETCH_OPTS \
static option_struct sketch_long_options[] = {\
    SHARED_OPTS\
};



namespace dashing2 {
void sketch_usage() {
    std::fprintf(stderr, "dashing2 sketch <opts> [fastas... (optional)]\n"
                         "We use only m-mers; if w <= k, however, this reduces to k-mers if the -w/--window-size is unspecified.\n"
                         SHARED_DOC_LINES
    );
}


int sketch_main(int argc, char **argv) {
    int c;
    int k = 16, w = -1, nt = -1;
    SketchSpace sketch_space = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = true;
    bool exact_kmer_dist = false, hpcompress = false;
    bool refine_exact = false;
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
    double downsample_frac = 1.;
    uint64_t seedseed = 0;
    Measure measure = SIMILARITY;
    std::ios_base::sync_with_stdio(false);
    std::string fsarg;
    // By default, use full hash values, but allow people to enable smaller
    bool normalize_bed = false;
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    std::string spacing;
    SKETCH_OPTS
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:CNs2BPWh?ZJGH", sketch_long_options, &option_index)) >= 0;) {
        switch(c) {
            SHARED_FIELDS
            case '?': case 'h': sketch_usage(); return 1;
        }
        //std::fprintf(stderr, "After getopt argument %d, of is %s\n",c , to_string(of).data());
    }
    const std::string ex(std::filesystem::absolute(std::filesystem::path(argv[-1])));
    std::string cmd(ex);
    for(char **s = argv; *s; cmd += std::string(" ") + *s++);
    std::fprintf(stderr, "[Dashing2] Invocation: %s ", cmd.data());
    if(nt < 0) {
        char *s = std::getenv("OMP_NUM_THREADS");
        if(s) nt = std::max(std::atoi(s), 1);
    }
    OMP_ONLY(omp_set_num_threads(nt));
    std::vector<std::string> paths(argv + optind, argv + argc);
    std::unique_ptr<std::vector<std::string>> qup;
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
    opts.downsample(downsample_frac);
    if(hpcompress) {
        if(!opts.homopolymer_compress_minimizers_) THROW_EXCEPTION(std::runtime_error("Failed to hpcompress minimizers"));
    }
    opts.filterset(fsarg);
    if((opts.sspace_ == SPACE_PSET || opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_EDIT_DISTANCE)
            && opts.kmer_result_ == ONE_PERM) {
        opts.kmer_result_ = FULL_SETSKETCH;
    }
    opts.bed_parse_normalize_intervals_ = normalize_bed;
    if(paths.empty()) {
        std::fprintf(stderr, "No paths provided. See usage.\n");
        sketch_usage();
        return 1;
    }
    SketchingResult result = sketch_core(opts, paths, outfile);
    result.nqueries(nq); // TODO: use nqueries to perform asymmetric comparisons
    if(cmpout.size()) {
        Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist, refine_exact);
        distopts.measure_ = measure;
        cmp_core(distopts, result);
    }
    return 0;
}

} // namespace dashing2
