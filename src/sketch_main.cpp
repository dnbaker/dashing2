#include "sketch_core.h"
#include "options.h"
#include "cmp_main.h"



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
    int k = -1, w = -1, nt = -1;
    SketchSpace sketch_space = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = true;
    bool exact_kmer_dist = false, hpcompress = false;
    bool refine_exact = false;
    long double compressed_a = -1.L, compressed_b = -1.L;
    bool fasta_dedup = false;
    double similarity_threshold = -1.;
    unsigned int count_threshold = 0.;
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
    int by_chrom = false;
    double downsample_frac = 1.;
    uint64_t seedseed = 0;
    size_t batch_size = 0;
    int nLSH = 2;
#if 0
    std::vector<std::pair<uint32_t, uint32_t>> compareids; // TODO: consider a sparse mode comparing only pairs of presented genomes.
#endif
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
            case OPTARG_HELP: case '?': case 'h': sketch_usage(); return 1;
        }
        //std::fprintf(stderr, "After getopt argument %d, of is %s\n",c , to_string(of).data());
    }
    if(k < 0) k = nregperitem(rht, use128);
    if(nt < 0) {
        char *s = std::getenv("OMP_NUM_THREADS");
        if(s) nt = std::max(std::atoi(s), 1);
    }
    OMP_ONLY(omp_set_num_threads(nt));
    std::vector<std::string> paths(argv + optind, argv + argc);
    std::unique_ptr<std::vector<std::string>> qup;
    if(ffile.size()) {
        if(!bns::isfile(ffile)) THROW_EXCEPTION(std::runtime_error("No path found at "s + ffile));
        std::ifstream ifs(ffile);
        static constexpr size_t bufsize = 1<<18;
        std::unique_ptr<char []> buf(new char[bufsize]);
        ifs.rdbuf()->pubsetbuf(buf.get(), bufsize);
        for(std::string l;std::getline(ifs, l);) {
            paths.push_back(l);
        }
        if(paths.empty()) {
            THROW_EXCEPTION(std::runtime_error("No paths read from "s + ffile));
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
        .count_threshold(count_threshold)
        .homopolymer_compress_minimizers(hpcompress)
        .seedseed(seedseed)
        .fasta_dedup(fasta_dedup);
    opts.by_chrom_ = by_chrom;
    opts.downsample(downsample_frac);
    opts.compressed_a_ = compressed_a;
    opts.compressed_b_ = compressed_b;
    opts.fd_level_ = nbytes_for_fastdists;
    opts.set_sketch_compressed();
    if(hpcompress) {
        if(!opts.homopolymer_compress_minimizers_) THROW_EXCEPTION(std::runtime_error("Failed to hpcompress minimizers"));
    }
    opts.filterset(fsarg);
    if((opts.sspace_ == SPACE_PSET || opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_EDIT_DISTANCE)
            && opts.kmer_result_ == ONE_PERM) {
        opts.kmer_result_ = FULL_SETSKETCH;
    }
    opts.bed_parse_normalize_intervals_ = normalize_bed;
    Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist, refine_exact, nLSH);
    if(paths.empty()) {
        std::fprintf(stderr, "No paths provided. See usage.\n");
        sketch_usage();
        return 1;
    }
    SketchingResult result;
    sketch_core(result, opts, paths, outfile);
    result.nqueries(nq);
    if(cmpout.size()) {
        distopts.measure_ = measure;
        distopts.cmp_batch_size_ = default_batchsize(batch_size, distopts);
        cmp_core(distopts, result);
    }
    return 0;
}

} // namespace dashing2
