#include "d2.h"

int main_usage() {
    std::fprintf(stderr, "dashing2 has several subcommands\n");
    std::fprintf(stderr, "Usage can be seen in those commands.\n");
    std::fprintf(stderr, "sketch: converts FastX into k-mer sets/sketches, and sketches BigWig and BED files\n");
    return 1;
}

int sketch_main(int argc, char **argv);

int main(int argc, char **argv) {
    if(argc == 1) {
        goto end;
    }
    if(std::strcmp(argv[1], "sketch") == 0)
        return sketch_main(argc - 1, argv + 1);
    end:
    return main_usage();
}
