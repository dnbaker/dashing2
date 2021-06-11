#include "fastxsketch.h"
#include "glob.h"


namespace dashing2 {

std::vector<std::string> glob(std::string pattern) {
    glob_t res;
    std::memset(&res, 0, sizeof(res));
    std::vector<std::string> ret;
    int rc = ::glob(pattern.data(), GLOB_TILDE, nullptr, &res);
    if(rc) throw std::runtime_error(std::string("glob error: ") + to_string(rc));
    ret.assign(res.gl_pathv, res.gl_pathv + res.gl_pathc);
    ::globfree(&res);
    return ret;
}
}

void usage() {
    std::fprintf(stderr, "readfq <opts> [fastas... (optional)]\n"
                         "We use only m-mers; if w <= k, however, this reduces to k-mers\n"
                         "-k: set k\n"
                         "-w: set window size for winnowing; by default, all m-mers are used.\n"
                         "-P: Sketch m-mers into ProbMinHash. Treats weighted sets as discrete probability distributions.\n"
                         "-B: Sketch m-mers into BagMinHash. Treats weighted sets as multisets.\n"
                         "-F: read paths from file in addition to positional arguments\n"
                         "-c: Use Count-Sketch counting instead of exact counting, using [arg] as the size.\n    "
                         "This accelerates weighted sketching at the cost of some approximation.\n"
                         "-s: Save m-mers. This puts the m-mers saved into .kmer files to correspond with the minhash samples.\n"
                         "-N: Save m-mer counts for sketches. This puts the m-mer counts saved into .kmercounts.f64 files to correspond with the m-mers.\n"
                         "-p: Set number of threads [1]\n"
                         "-Z: Full setsketch (not stochastically-averaged)\n"
                         "-H: Full m-mer set. This generates a sorted hash set for m-mers in the data. If the parser is windowed (-w is set), this may be rather small.\n"
                         "-J: Full m-mer countdict. This generates a sorted hash set for m-mers in the data, and additionally saves the associated counts for these m-mers.\n"
                         "-G: Full m-mer sequence.\n"
                         "-2: Use 128-bit k-mer hashes instead of 64-bit\n"
    );
}
using namespace dashing2;
int main(int argc, char **argv) {
    int c;
    int k = 16, w = 50, nt = 1;
    SketchSpace s = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false;
    size_t cssize = 0, sketchsize = 1024;
    std::string ffile;
    for(;(c = getopt(argc, argv, "p:k:w:c:f:S:Ns2SBPWh?ZJGH")) >= 0;) {switch(c) {
        case 'k': k = std::atoi(optarg); break;
        case 'w': w = std::atoi(optarg); break;
        case 'W': cache = true; break;
        case 'B': s = SPACE_MULTISET; res = FULL_SETSKETCH; break;
        case 'P': s = SPACE_PSET; res = FULL_SETSKETCH; break;
        case 'Z': res = FULL_SETSKETCH; break;
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;
        case 'p': nt = std::atoi(optarg); break;
        case 'S': sketchsize = std::atoi(optarg); break;
        case 'N': save_kmercounts = true; [[fallthrough]];
        case 's': save_kmers = true; break;
        case 'H': res = FULL_MMER_SET; break;
        case 'J': res = FULL_MMER_COUNTDICT; break;
        case 'G': res = FULL_MMER_SEQUENCE; break;
        case 'F': ffile = optarg; break;
        case '2': use128 = true; break;
        case '?': case 'h': usage(); std::exit(1);
    }}
    std::vector<std::string> paths(argv + optind, argv + argc);
    if(ffile.size()) {
        std::ifstream ifs(ffile);
        for(std::string l;std::getline(ifs, l);) {
            paths.push_back(l);
        }
    }
    Dashing2Options opts(k, w, bns::DNA, s);
    opts
        .nthreads(nt)
        .kmer_result(res)
        .cache_sketches(cache)
        .cssize(cssize)
        .use128(use128)
        .sketchsize(sketchsize);
    opts.save_kmers_ = save_kmers;
    opts.save_kmercounts_ = save_kmercounts;
    if(paths.empty()) {
        paths = glob("bonsai/test/*.fna.gz");
        std::fprintf(stderr, "Instead, sketch the bonsai/test/*fna.gz test files\n");
    }
    FastxSketchingResult result = fastx2sketch(opts, paths);
    result.print();
}
