#include "bwsketch.h"
using namespace dashing2;

#define gett std::chrono::high_resolution_clock::now

int main(int argc, char **argv) {
    if(argc != 2) {
        std::fprintf(stderr, "Usage: readbw <in.bw>\n");
        return 1;
    }
    Dashing2Options opts(1, -1, bns::DNA, SPACE_PSET, BIGWIG);
    int nt = 16;
    if(char *s = std::getenv("OMP_NUM_THREADS")) nt = std::atoi(s);
    opts.nthreads(nt);
    auto t = gett();
    auto s = bw2sketch(argv[1], opts, /*parallel_process=*/true);
    auto tend = gett();
    std::fprintf(stderr, "Sketching bigwig file took %gms\n", std::chrono::duration<double, std::milli>(tend - t).count());
    for(const auto &pair: *s.chrmap_) {
        std::fprintf(stderr, "key %s\n", pair.first.data());
    }
    return 0;
}
