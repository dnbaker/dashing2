#include "src/d2.h"
#include "minispan.h"
#include "mio.hpp"
#include "fmt/format.h"

namespace dashing2 {
int printmin_main(int argc, char **argv) {
    bool emit_fasta = false;
    std::FILE *ofp = stdout;
    if(auto it = std::find_if(argv, argv + argc, [](auto x) {return !(std::strcmp(x, "-h") && std::strcmp(x, "--help"));}); it != argv + argc) {
        goto usage;
    }
    for(int c;(c = getopt(argc, argv, "fo:h")) >=0;) {switch(c) {
        case 'f':
            emit_fasta = true; break;
        case 'o': ofp = std::fopen(optarg, "w"); if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open ") + optarg));
                break;
        case 'h':
            usage:
            std::fprintf(stderr, "Usage: -f: emit fasta. Default - emits tabular result.\n-o: Write to file instead of stdout.\n");
            return 1;
    }}
    buffer_to_blksize(ofp);
    if(optind == argc) THROW_EXCEPTION(std::invalid_argument("Required: one positional argument for printmin_main"));
    std::string inpath = argv[optind];
    mio::mmap_source ifp(inpath);
    size_t nseqs;
    std::memcpy(&nseqs, ifp.data(), sizeof(nseqs));
    uint32_t k, w, dtype;
    {
        auto ptr = (uint32_t *)(ifp.data() + 8);
        k = ptr[0];
        w = ptr[1];
        dtype = ptr[2];
    }
    //const bool canon = (dtype >> 8) & 1;
    const bns::RollingHashingType rht = static_cast<bns::RollingHashingType>(dtype & 0xff);
    if(rht != bns::InputType::DNA) {
        THROW_EXCEPTION(std::runtime_error("Not yet implemented: minimizer sequence printing for non-DNA alphabets"));
    }
    const double *lptr = (const double *)(ifp.data() + sizeof(nseqs) + sizeof(uint32_t) * 3);
    const minispan<double> lspan(lptr, nseqs);
    std::vector<size_t> offsets(nseqs + 1);
    for(size_t i = 0; i < nseqs; ++i) {
        offsets[i + 1] = offsets[i] + lspan[i];
    }
    const size_t totalsum = offsets.back();
    size_t totalfilesize = sizeof(size_t)
                           + sizeof(uint32_t) * 3
                           + nseqs * sizeof(double)
                           + totalsum * sizeof(uint64_t);
    if(totalfilesize != ifp.size()) {
        THROW_EXCEPTION(std::runtime_error(std::string("Unexpected filesize ") + std::to_string(ifp.size()) + " vs " + std::to_string(totalfilesize) + ". Corrupted file?"));
    }
    const uint64_t *kmerptr = (const uint64_t *)(lptr + nseqs);
    const minispan<uint64_t> kspan(kmerptr, totalsum);
    bns::Spacer sp(k, w);
    for(size_t seqid = 0; seqid < nseqs; ++seqid) {
        const minispan<uint64_t> subminispan(&kspan[offsets[seqid]], offsets[seqid + 1] - offsets[seqid]);
        if(emit_fasta) {
            for(size_t i = 0; i < subminispan.size(); ++i) {
                fmt::print(ofp, ">MinimizerSequence{}-Minimizer#{}\n{}\n", seqid, i, sp.to_string(subminispan[i]));
            }
        } else {
            fmt::print(ofp, "MinimizerSequence{}", seqid);
            for(const auto v: subminispan) {
                fmt::print(ofp, " {}", sp.to_string(v));
            }
            fmt::print(ofp, "\n");
        }
    }
#if 0
    for(const int32_t L: lspan) {
        if(emit_fasta) {
            size_t minid = 0;
            for(size_t i = 0; i < L; ++i) {
                uint64_t item;
                checked_fread(&item, 1, sizeof(item), ifp);
                fmt::print(ofp, ">MinimizerSequence{}-Minimizer#{}\n{}\n", id, minid++, sp.to_string(item));
            }
        } else {
            fmt::print(ofp, "MinimizerSequence{}", id);
            for(size_t i = 0; i < L; ++i) {
                uint64_t item;
                checked_fread(&item, 1, sizeof(item), ifp);
                fmt::print(ofp, " {}", sp.to_string(item));
            }
            fmt::print(ofp, "\n");
        }
        ++id;
    }
#endif
    return 0;
}

}
