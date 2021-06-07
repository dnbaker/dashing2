#include "d2.h"

namespace dashing2 {
int cmp_main(int argc, char **argv);
int sketch_main(int argc, char **argv);
} // dashing2

int main_usage() {
    std::fprintf(stderr, "dashing2 has several subcommands\n");
    std::fprintf(stderr, "Usage can be seen in those commands.\n");
    std::fprintf(stderr, "cmp: compres previously sketched/decomposed k-mer sets and emits results.\n");
    std::fprintf(stderr, "sketch: converts FastX into k-mer sets/sketches, and sketches BigWig and BED files; also contains functionality from cmp, for one-step sketch and comparisons\n");
    return 1;
}
using namespace dashing2;


int main(int argc, char **argv) {
    if(argc > 1) {
        if(std::strcmp(argv[1], "sketch") == 0)
            return sketch_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "cmp") == 0 || std::strcmp(argv[1], "dist") == 0)
            return cmp_main(argc - 1, argv + 1);
    }
    return main_usage();
}
