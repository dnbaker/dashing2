#include "fastxsketch.h"

#if 0
    ParseOptions(int k, int w=-1, bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, CountingType count=EXACT_COUNTING, DataType dtype=FASTX, size_t nt=0, bool use128=false):
enum KmerSketchResultType {
    ONE_PERM = 0,       // Faster (3-4x) than Full, comparable accuracy for both cardinality and set similarities
    // This is a stochastically-averaged generalized HyperLogLog
    // Constant-time updates,
    FULL_SETSKETCH = 1, // Not stochastically-averaged; potentially better LSH properties
    // This is a generalized HyperLogLog
    FULL_MMER_SET = 2,
    /*
     * Convert the genome into a k-mer set; uses a hash table
    */
    FULL_MMER_SEQUENCE = 3,
    /*
    Convert the genome into a list of k-mers/minimizers; could be used for minimizer index generation
    */
    FULL_MMER_COUNTDICT = 4
    /*
        Convert into a k-mer: count dictionary.
    */
};

#endif
using namespace dashing2;
int main(int argc, char **argv) {
    int c;
    int k = 16, w = 50, nt = 1;
    SketchSpace s = SPACE_SET;
    CountingType ct = EXACT_COUNTING;
    KmerSketchResultType res = ONE_PERM;
    for(;(c = getopt(argc, argv, "p:k:w:BPC")) >= 0;) {switch(c) {
        case 'k': k = std::atoi(optarg); break;
        case 'w': w = std::atoi(optarg); break;
        case 'B': s = SPACE_MULTISET; res = FULL_SETSKETCH; break;
        case 'P': s = SPACE_PSET; res = FULL_SETSKETCH; break;
        case 'c': ct = COUNTSKETCH_COUNTING; break;
        case 'p': nt = std::atoi(optarg); break;
    }}
    ParseOptions opts(k, w, bns::DNA, s, ct);
    opts.nthreads(nt);
    opts.kmer_result(res);
    std::vector<std::string> paths(argv + optind, argv + argc);
    if(paths.empty()) {
        std::fprintf(stderr, "Usage: readfx (flags) <path> <path2> ...\n");
        std::exit(1);
    }
    FastxSketchingResult result = fastx2sketch(opts, paths);
}
