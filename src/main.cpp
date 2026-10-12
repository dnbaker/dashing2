#include "d2.h"
#include <cstdio>
#include <cstring>
#include <filesystem>

namespace dashing2 {
int cmp_main(int argc, char **argv);
int contain_main(int argc, char **argv);
int wsketch_main(int argc, char **argv);
int sketch_main(int argc, char **argv);
int printmin_main(int argc, char **argv);
}

int main_usage() {
    std::fprintf(stderr, "dashing2 has several subcommands: sketch, cmp, wsketch, and contain.\n");
    std::fprintf(stderr, "Usage can be seen in those subcommands. (e.g., `dashing2 sketch -h`)\n\n");
    std::fprintf(stderr, "\tsketch: converts FastX into k-mer sets/sketches, and sketches BigWig and BED files; also contains functionality from cmp, for one-step sketch and comparisons\n"
                         "This is probably the most common subcommand to use.\n\n"
    );
    std::fprintf(stderr, "\tcmp: compares previously sketched/decomposed k-mer sets and emits results. alias: dist\n\n");
    std::fprintf(stderr, "\tcontain: Takes a k-mer database (built with dashing2 sketch --save-kmers), then computes coverage for all k-mer references using input streams.\n");
    std::fprintf(stderr, "\twsketch: Takes a tuple of [1-3] input binary files [(u32 or u64), (float or double), (u32 or u64)] and performs weighted minhash sketching.\n"
                         "Three files are treated as Compressed Sparse Row (CSR)-format, where the third file contains indptr values, specifying the lengths of consecutive runs of pairs in the first two files corresponding to each row.\n"
                         "wsketch is for sketching binary files which have already been summed, whereas sketch is for parsing and sketching (from Fast{qa}, BED, BigWig)\n");
    std::fprintf(stderr, "\n\nMiscellania:\n");
    std::fprintf(stderr, "printmin: Emit minimizer sequence sets in human-readable form.\n");
    return 1;
}
using namespace dashing2;


int main(int argc, char **argv) {
    std::string cmd(std::filesystem::absolute(std::filesystem::path(argv[0])));
    for(char **s = (argv + 1); *s; cmd += std::string(" ") + *s++);
    std::fprintf(stderr, "#Calling Dashing2 version %s with command '%s'\n", DASHING2_VERSION, cmd.data());
    if(argc > 1) {
        if(std::strcmp(argv[1], "sketch") == 0)
            return sketch_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "cmp") == 0 || std::strcmp(argv[1], "dist") == 0)
            return cmp_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "wsketch") == 0)
            return wsketch_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "contain") == 0)
            return contain_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "printmin") == 0) {
            return printmin_main(argc - 1, argv + 1);
        }
    }
    return main_usage();
}
