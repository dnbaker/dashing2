#include "bedsketch.h"
using namespace dashing2;

#define gett std::chrono::high_resolution_clock::now

int main(int argc, char **argv) {
    if(argc != 2) {
        std::fprintf(stderr, "Usage: readbed <in.bed>\n");
        return 1;
    }
    Dashing2Options opts(1, -1, bns::DNA, SPACE_PSET, BED);
    int nt = 16;
    if(char *s = std::getenv("OMP_NUM_THREADS")) nt = std::atoi(s);
    opts.nthreads(nt);
    opts.sketchsize(64);
    auto t = gett();
    auto s = bed2sketch(argv[1], opts);
    auto tend = gett();
    std::fprintf(stderr, "Sketching BED file took %gms\n", std::chrono::duration<double, std::milli>(t - tend).count());
    opts.cssize(10000);
    t = gett();
    s = bed2sketch(argv[1], opts);
    tend = gett();
    std::fprintf(stderr, "Sketching BED file with count sketch took %gms\n", std::chrono::duration<double, std::milli>(t - tend).count());
    opts.cssize(0);
    opts.sspace_ = SPACE_SET;
    t = gett();
    s = bed2sketch(argv[1], opts);
    tend = gett();
    std::fprintf(stderr, "Sketching BED file in set space took %gms\n", std::chrono::duration<double, std::milli>(t - tend).count());
    return 0;
}
