#pragma once
#ifndef DASHING2_OPTIONS_H__
#define DASHING2_OPTIONS_H__
#include <enums.h>
#include <getopt.h>

namespace dashing2 {
#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

using option_struct = struct option;

enum OptArg {
    OPTARG1 = 1000,
    OPTARG_BED,
    OPTARG_BIGWIG,
    OPTARG_LEAFCUTTER,
    OPTARG_OUTPREF,
    OPTARG_CMPOUT,
    OPTARG_BINARY_OUTPUT,
    OPTARG_FASTCMP,
    OPTARG_REGBYTES,
    OPTARG_BBIT_SIGS,
    OPTARG_EXACT_KMER_DIST,
    OPTARG_ISZ,
    OPTARG_BED_NORMALIZE,
    OPTARG_PROTEIN,
    OPTARG_PROTEIN6,
    OPTARG_PROTEIN8,
    OPTARG_PROTEIN14,
    OPTARG_HPCOMPRESS,
    OPTARG_DOWNSAMPLE_FRACTION,
    OPTARG_REFINEEXACT,
    OPTARG_MASHDIST,
    OPTARG_SYMCONTAIN,
    OPTARG_CONTAIN,
    OPTARG_ASYMMETRIC_ALLPAIRS,
    OPTARG_SET,
    OPTARG_SPACING,
    OPTARG_RANDOM_SEED,
    OPTARG_FILTERSET,
    OPTARG_PARSEBYSEQ,
    OPTARG_HELP,
    OPTARG_CMP_BATCH_SIZE,
    OPTARG_GREEDY,
    OPTARG_NLSH,
    OPTARG_ENTROPYMIN,
    OPTARG_SIGRAMLIMIT,
    OPTARG_DUMMY
};

#define SHARED_OPTS \
    LO_ARG("ffile", 'F')\
    LO_ARG("qfile", 'Q')\
    LO_ARG("threads", 'p')\
    LO_ARG("sketchsize", 'S')\
    LO_ARG("cmpout", OPTARG_CMPOUT)\
    LO_ARG("distout", OPTARG_CMPOUT)\
    LO_ARG("cmp-outfile", OPTARG_CMPOUT)\
    LO_ARG("outprefix", OPTARG_OUTPREF)\
    LO_ARG("topk", 'K')\
    LO_ARG("similarity-threshold", 'T')\
    LO_ARG("fastcmp", OPTARG_FASTCMP)\
    LO_ARG("countsketch-size", 'c')\
    LO_ARG("window-size", 'w')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("outfile", 'o')\
    LO_ARG("count-threshold", 'm')\
    LO_ARG("threshold", 'm')\
    LO_FLAG("binary-output", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("probs", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("pminhash", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("pmh", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("PMH", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("probminhash", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
    /*LO_FLAG("edit-distance", 'E', sketch_space, SPACE_EDIT_DISTANCE)*/\
    /*LO_FLAG("set", 'H', res, FULL_MMER_SET)*/\
    /*LO_FLAG("exact-kmer-dist", OPTARG_EXACT_KMER_DIST, exact_kmer_dist, true)*/\
    LO_FLAG("bbit-sigs", OPTARG_BBIT_SIGS, truncate_mode, 1)\
    LO_FLAG("intersection", OPTARG_ISZ, measure, INTERSECTION)\
    LO_FLAG("mash-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("symmetric-containment", OPTARG_SYMCONTAIN, measure, SYMMETRIC_CONTAINMENT)\
    LO_FLAG("containment", OPTARG_CONTAIN, measure, SYMMETRIC_CONTAINMENT)\
    LO_FLAG("poisson-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("compute-edit-distance", OPTARG_MASHDIST, measure, M_EDIT_DISTANCE)\
    LO_FLAG("multiset", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("bmh", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("BMH", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
    LO_FLAG("asymmetric-all-pairs", OPTARG_ASYMMETRIC_ALLPAIRS, ok, OutputKind::ASYMMETRIC_ALL_PAIRS)\
    LO_FLAG("phylip", OPTARG_PHYLIP, ok, OutputKind::PHYLIP)\
    LO_ARG("regbytes", OPTARG_REGBYTES)\
    /*LO_ARG("set", 'H')*/\
    {"save-kmers", no_argument, 0, 's'},\
    {"save-kmercounts", no_argument, 0, 'N'},\
    {"hp-compress", no_argument, 0, OPTARG_HPCOMPRESS},\
    {"refine-exact", no_argument, 0, OPTARG_REFINEEXACT},\
    {"edit-distance", no_argument, 0, 'E'},\
    {"oneperm-setsketch", no_argument, 0, 'Z'},\
    {"oneperm", no_argument, 0, 'Z'},\
    {"one-perm", no_argument, 0, 'Z'},\
    {"oph", no_argument, 0, 'Z'},\
    {"doph", no_argument, 0, 'Z'},\
    {"normalize-intervals", no_argument, 0, OPTARG_BED_NORMALIZE},\
    {"protein", no_argument, 0, OPTARG_PROTEIN},\
    {"protein20", no_argument, 0, OPTARG_PROTEIN},\
    {"enable-protein", no_argument, 0, OPTARG_PROTEIN},\
    {"protein6", no_argument, 0, OPTARG_PROTEIN6},\
    {"protein8", no_argument, 0, OPTARG_PROTEIN8},\
    {"protein14", no_argument, 0, OPTARG_PROTEIN14},\
    {"downsample", required_argument, 0, OPTARG_DOWNSAMPLE_FRACTION},\
    {"cache", no_argument, 0, 'W'},\
    {"no-canon", no_argument, 0, 'C'},\
    {"set", no_argument, 0, OPTARG_SET},\
    {"exact-kmer-dist", no_argument, 0, OPTARG_EXACT_KMER_DIST},\
    {"spacing", required_argument, 0, OPTARG_SPACING},\
    {"seed", required_argument, 0, OPTARG_RANDOM_SEED},\
    {"filterset", required_argument, 0, OPTARG_FILTERSET},\
    {"parse-by-seq", no_argument, 0, OPTARG_PARSEBYSEQ},\
    {"help", no_argument, 0, OPTARG_HELP},\
    {"batch-size", required_argument, 0, OPTARG_CMP_BATCH_SIZE},\
    {"greedy", required_argument, 0, OPTARG_GREEDY},\
    {"nlsh", required_argument, 0, OPTARG_NLSH},\
    {"nLSH", required_argument, 0, OPTARG_NLSH},\
    {"entmin", no_argument, 0, OPTARG_ENTROPYMIN},\
    {"by-chrom", no_argument, (int *)&by_chrom, 1},\
    {"sketch-size-l2", required_argument, 0, 'L'},\
    {"sig-ram-limit", required_argument, 0, OPTARG_SIGRAMLIMIT}



#define TOPK_FIELD case 'K': {ok = OutputKind::KNN_GRAPH; topk_threshold = std::atoi(optarg); break;}
#define SIMTHRESH_FIELD case 'T': {ok = OutputKind::NN_GRAPH_THRESHOLD; similarity_threshold = std::atof(optarg); break;}
#define GREEDY_FIELD case OPTARG_GREEDY: {ok = OutputKind::DEDUP; similarity_threshold = std::atof(optarg); break;}
#define CMPOUT_FIELD case OPTARG_CMPOUT: {cmpout = optarg; break;}
#define FASTCMP_FIELD case OPTARG_FASTCMP: {nbytes_for_fastdists = std::atof(optarg);\
            if(nbytes_for_fastdists != 8. && nbytes_for_fastdists != 4. && nbytes_for_fastdists != 2. && nbytes_for_fastdists != 1. && nbytes_for_fastdists != .5){\
                std::fprintf(stderr, "--fastcmp must have 8, 4, 2, 1, or 0.5 as the argument. These are the only register sizes supported.\n");\
                throw std::runtime_error("See usage for --fastcmp instructions.");\
            }\
            break;\
            }
#define PROT_FIELD case OPTARG_PROTEIN: {rht = bns::PROTEIN20; canon = false; std::fprintf(stderr, "Parsing 20-character amino acod sequences\n"); break;} \
    case OPTARG_PROTEIN6: {rht = bns::PROTEIN_6; canon = false; std::fprintf(stderr, "Parsing amino acid sequences with 6-letter compressed alphabet.\n"); break;}\
    case OPTARG_PROTEIN14: {rht = bns::PROTEIN14; canon = false; std::fprintf(stderr, "Parsing amino acid sequences with 14-letter compressed alphabet.\n"); break;}\
    case OPTARG_PROTEIN8: {rht = bns::PROTEIN8; canon = false; std::fprintf(stderr, "Using 3-bit protein encoding\n"); break;}
#define REFINEEXACT_FIELD case OPTARG_REFINEEXACT: {refine_exact = true; break;}

#define SHARED_FIELDS TOPK_FIELD SIMTHRESH_FIELD \
        CMPOUT_FIELD FASTCMP_FIELD PROT_FIELD REFINEEXACT_FIELD \
        GREEDY_FIELD \
        case 'E': sketch_space = SPACE_EDIT_DISTANCE; break;\
        case 'C': canon = false; break;\
        case 'p': nt = std::atoi(optarg); break;\
        case OPTARG_EXACT_KMER_DIST: {\
            exact_kmer_dist = true;\
            std::fprintf(stderr, "Exact kmer distance is set to %d -- %d/%d\n", exact_kmer_dist, c, option_index);\
            break;\
        }\
        case 'S': sketchsize = std::atoi(optarg); break;\
        case 'N': save_kmers = save_kmercounts = true; break;\
        case 's': save_kmers = true; break;\
        case OPTARG_SET: case 'H': res = FULL_MMER_SET; break;\
        case 'J': res = FULL_MMER_COUNTDICT; break;\
        case 'G': res = FULL_MMER_SEQUENCE; break;\
        case '2': use128 = true; break;\
        case 'm': count_threshold = std::atoi(optarg); break;\
        case 'F': ffile = optarg; break;\
        case 'Q': qfile = optarg; ok = PANEL; break;\
        case OPTARG_BED_NORMALIZE: normalize_bed = true; break;\
        case OPTARG_REGBYTES: nbytes_for_fastdists = std::atof(optarg); break;\
        case 'o': outfile = optarg; break;\
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;\
        case 'k': k = std::atoi(optarg); break;\
        case 'w': w = std::atoi(optarg); break;\
        case 'W': cache = true; break;\
        case 'B': sketch_space = SPACE_MULTISET; res = FULL_SETSKETCH; break;\
        case 'P': sketch_space = SPACE_PSET; res = FULL_SETSKETCH; break;\
        case 'Z': res = ONE_PERM; break;\
        case OPTARG_ISZ: measure = INTERSECTION; break;\
        case OPTARG_OUTPREF: {\
            outprefix = optarg; break;\
        } \
        case OPTARG_HPCOMPRESS: {\
            hpcompress = true; break;\
        }\
        case OPTARG_DOWNSAMPLE_FRACTION: {\
            downsample_frac = std::atof(optarg); break;\
        }\
        case OPTARG_SPACING: spacing = optarg; canon = false; break;\
        case OPTARG_RANDOM_SEED: {seedseed = std::strtoull(optarg, 0, 10);} break;\
        case OPTARG_FILTERSET: fsarg = optarg; break;\
        case OPTARG_PARSEBYSEQ: parse_by_seq = true; break;\
        case OPTARG_CMP_BATCH_SIZE: batch_size = std::strtoull(optarg, 0, 10); break;\
        case OPTARG_NLSH: nLSH = std::atoi(optarg); break;\
        case OPTARG_ENTROPYMIN: entmin = true; break;\
        case 'L': {\
            const int ssl2 = std::atoi(optarg);\
            sketchsize = size_t(1) << ssl2;\
            if(ssl2 <= 0 || ssl2 >= 64) {\
                std::fprintf(stderr, "Error: ssl2 is out of bounds. ssl2: %d. sketchsize: %zu. Did you mean to specify sketchsize \n", ssl2, sketchsize);\
                std::exit(1);\
            }\
            std::fprintf(stderr, "Using log_2 sketchsize = %d, yielding sketchsize = %zu\n", ssl2, sketchsize);\
            break;\
        }\
        case OPTARG_SIGRAMLIMIT: {\
            MEMSIGTHRESH = std::strtoull(optarg, nullptr, 10);\
        } break;



static constexpr const char *siglen =
        sizeof(RegT) == 2 ? "16":
        sizeof(RegT) == 1 ? "8":
        sizeof(RegT) == 4 ? "32":
        sizeof(RegT) == 8 ? "64": "128";
#define SHARED_DOC_LINES \
        "Flags:\n\n"\
        "Runtime options:\n"\
        "-p/--threads: Set number of threads [1]\n"\
        "Encoding options\n"\
        "Dashing2 can sketch 4 kinds of files:"\
        "Fastq/Fasta, which has specific encoding options (default)\n"\
        "--bed to sketch BED files for interval sets\n"\
        "--bigwig to sketch BigWig files for coverage vectors\n"\
        "and --leafcutter to sketch LeafCutter splicing output\n"\
        "\n\nFastx Options:\n"\
        "-k/--kmer-length: set k. Defaults to the largest k expressible directly in uint64_t.\n"\
        "If k is greater than this limit (31 for DNA, 14 for --protein, 22 for --protein8, 24 for --protein6), then rolling hashes will be generated instead of exact k-mer encodings.\n"\
        "-w/--window-size: set window size for winnowing; by default, all k-mers are used. If w > k, then only the minimum-hash k-mer in each window is processed\n"\
        "This can be useful to speed up sketching or to reduce the number of items sketched in exact mode.\n"\
        "--entmin: If -w/--window-size is enabled, this option weights the hash value by the entropy of the k-mer itself.\nThis is only valid for k-mers short enough to be encoded exactly in 64-bit or 128-bit integers, depending on if --long-kmers is enabled.\n"\
        "--spacing: Set a spacing scheme for spaced minimizers\n"\
        "This must have 1 less integer than the k-mer length.\n"\
        "e.g., -k5 --spacing 0,1,1,0 specifies a match length of 5 with a match pattern of `KK$K$KK`, where $ positions are ignored and `K` positions are kept.\n"\
        "This can also be run-length compressed in <space, num> format.\n"\
        "For example, --spacing 0,1x2,0 is equivalent to 0,1,1,0.\n"\
        "-2/--128bit/long-kmers: Use 128-bit k-mer hashes instead of 64-bit\n"\
        "-m/--threshold/--count-threshold <arg>: Set a count threshold for inclusion. Default: 0. If > 1, this will only sketch k-mers with count >= <arg>\n"\
        "\nFastx Alphabet:\n"\
        "Dashing2 sketches DNA be default. This can be changed with the following flags; this will disable canonicalization.\n"\
        "--enable-protein: Use 20 character amino acid alphabet.\n"\
        "--protein14: Use 14 character amino acid alphabet.\n"\
        "--protein6: Use 6 character amino acid alphabet.\n"\
        "--protein8: Use 8 character (3-bit) amino acid alphabet.\n"\
        "--no-canon: If DNA is being encoded, this disables canonicalization. By default, DNA sequence is canonicalized with its reverse-complement.\n"\
        "            Otherwise, this is ignored\n"\
        "--seed: Set a seed for k-mer hashing; If 0, this disables k-mer XORing and k-mers are encoded directly if a k-mer type can represent it.\n"\
        "        Otherwise, this changes the hash function applied to k-mers when generated sorted hash sets. This makes it easy to decode quickly, but we can still get good bottom-k estimates using these hashes\n"\
        "        the xor value for u64 kmers (unless --long-kmers is enabled) is the Wang 64-bit hash of the seed.\n"\
        "        u128 kmers (--long-kmers) have the same lower 64 bits, but the upper 64 bits are the Wang 64-bit hash of the u64 xor value.\n"\
        "\nPathsOptions\n\n"\
        "By default, dashing2 reads positional arguments and sketches them. You may want to use flags instructing it\n"\
        "to read from paths in <file>. Additionally, you can put multiple files separated by spaces into a single line "\
        "to place them all into a single sketch instead.\n"\
        "-F/--ffile: read paths from file in addition to positional arguments.\n"\
        "-Q/--qfile: read query paths from file; this is used for asymmetric queries (e.g., containment).\n"\
        "If multiple files are space-delimited in a single line, they will be sketched jointly.\n"\
        "K-mer Filtering\n\n"\
        "If there are a set of common k-mers or artefactual sequence, you can specify --filterset to skip k-mers in this file when sketching other files.\n"\
        "By default, this converts it into a sorted hash set and skips k-mers which are found in the set.\n"\
        "`--filterset [path]` yields this.\n"\
        "--downsample\t Downsample minimizers at fraction <arg> . Default is 1: IE, all minimizers pass.\n"\
        "\nSketch options\n"\
        "These decide how k-mers are accumulated.\n"\
        "Default behavior is set sketching (tossing multiplicities). If --multiset or --prob is set or a minimum count is provided,"\
        "\nk-mers will be counted before sketching.\n"\
        "-S/--sketchsize: Set sketchsize (1024)\n"\
        "Warning: this is a significant change from Dashing1, where sketch size was required to be a power of 2, and so the sketch size would be 2**<arg>.\n"\
        "This restriction is no longer in place, and so any positive integer can be selected.\n"\
        "e.g. - instead of -S12 in Dashing, you would specify -S 4096\n"\
        "For convenience, we offer another argument to specify the size in log2.\n"\
        "-L/--sketch-size-l2: Set sketchsize to 2^<arg>. Must be > 0 and < 64\n"\
        "\nSketching modes\n\n"\
        "When sketching, there are several sketch options.\n"\
        "We use the SetSketch as the default sketch structure; This ignores multiplicities, except for filtering \n"\
        "SetSketching -\n"\
        "1. FullSetSketch (the default), sketches in set space. Multiplicities are not considered.\n"\
        "This uses increasing exponential random variates, and is more reliable but slower than the One-Permutation.\n"\
        "2. One-Permutation (stochastically-averaged) SetSketch.\n"\
        "-Z/--doph/--oneperm/--oneperm-setsketch: \n"\
        "Stochastically-averaged setsketch. This is faster at sketching, but has a small probability of failure which grows with sketch size. Big one-permutation sketches may perform poorly.\n"\
        "This is faster to sketch, equal speed to compare; it may fail if the sketch size is too large;\n"\
        "If sketchsize is significantly smaller than the entities being processed, this should work;\n"\
        "This should perform similarly to default setsketch behavior, but has better behaviors with large sketches and small sets\n"\
        "Weighted Sketching -\n"\
        "We provide to weighted sketching algorithms -- BagMinHash and ProbMinHash\n"\
        "BagMinHash is an LSH for the weighted Jaccard similarity, which treats k-mer counts as weighted sets.\n"\
        "ProbMinHash is an LSH for the probability Jaccard index, which normalizes observations by total counts, yielding a discrete probability distribution for each collection.\n"\
        "Both require counting for sequence data, but do not for input methods with counting already performed, e.g. BigWig.\n"\
        "ProbMinHash is most applicable for datasets where sampling fractions are important, such as expression or splicing counts.\n"\
        "--prob: Sketch k-mers into ProbMinHash. Treats weighted sets as discrete probability distributions.\n"\
        "        Aliases: --pminhash, --probs, --pmh, --PMH\n"\
        "BagMinHash is most applicable for datasets where the absolute quantities of an event matter, rather than its fraction of the complete dataset.\n"\
        "-B/--multiset: Treats weighted sets as multisets.\n"\
        "        Aliases: --bagminhash, --bmh, --BMH\n"\
        "It typically comes at 2-4x runtime cost over ProbMinHash, depending on sketch size.\n"\
        "Weighted Sketching (Counting) Parameters -\n"\
        "-c/--countsketch-size: Use Count-Sketch counting instead of exact counting, using <arg> as the size.\n    "\
        "This allows you to avoid unbounded dictionary size at the cost of some approximation of the weighted sets\n"\
        "(This only affects BagMinHash and ProbMinHash, and only for sequence data.\n"\
        "\n\nExact Modes\n"\
        "You can also emit full k-mer sets, a count dictionary (key-count map), or a set of minimizer sequences\n"\
        "-H/--set: Full k-mer set. This generates a sorted hash set for k-mers in the data. If the parser is windowed (-w is fairly large), this could even be rather small.\n"\
        "                If an LSH table is generated, then weighted bottom-k hashes are used to build an LSH table\n"\
        "-J/--countdict: Full k-mer countdict. \n"\
        "                This generates a sorted hash set for k-mers in the data, and additionally saves the associated counts for these k-mers.\n"\
        "                If an LSH table is generated, then weighted bottom-k hashes as in Cohen, E. \"Summarizing Data using Bottom-K Sketches\"\n"\
        "-G/--seq: Full k-mer (or minimizer) sequence. This faster than building the hash set, and can be used to build a minimizer index afterwards\n"\
        "          On the other hand, it can require higher memory for large sequence collections\n"\
        "          If you use --parse-by-seq with this and an output path is provided, then the stacked minimizer sequences will be written to it.\n"\
        "          The format is the similar to the standard stacked sketches, except that the cardinality fields instead represent minimizer sequence lengths (in 64-bit registers).\n"\
        "          Specifically, it consists of a header: [uint64_t nitems, uint32_t k, uint32_t w], followed by `nitems` [double], specifying sequence lengths of 64-bit registers\n"\
        "    Dependent option (only for --seq/-G parsing)\n"\
        "          --hp-compress:\n"\
        "              Minimizer sequence will be homopolymer-compressed before emission. \n"\
        "              This makes the sequences insensitive to the lengths of minimizer stretches, which may simplify match finding\n"\
        "\n\nDistance Options (shared)\n"\
        "--Exhaustive Distance Outputs--\n"\
        "The default output format is all-pairs upper-triangular.\n"\
        "We provide three other full matrix methods - \n"\
        "1. PHYLIP Upper Triangular distance matrix, enabled by --phylip.\n"\
        "2. Square distance matrix, enabled by --asymmetric-all-pairs\n"\
        "--asymmetric-all-pairs: emit square distance matrix\n"\
        "Emits a full square distance matrix rather than upper-triangular.\n"\
        "3. Rectangular matrix, enabled by -Q/-qfile\n"\
        "Use this if you want to compare a set of queries to a set of references rather than complete all-pairs. Note: -F must be provided, or reference files should be added as positional arguments\n"\
        "positional arguments and -F paths are treated as a reference set;\n"\
        "Paths provided in -Q/--qfile are are treated as a query set.\n"\
        "The output shape then has |F| rows and |Q| columns, (F, Q) in row-major format\n"\
        "--Sparse Distance Outputs--\n"\
        "We provide several sparse distance options - \n"\
        "    1. topk, -- only the top-k nearest neighbors are emitted for each item\n"\
        "    2. thresholded, -- only the similarities > threshold are emitted\n"\
        "Both change the output format from a full matrix into compressed-sparse row (CSR) format listing only the filtered entries;\n"\
        "All of these are powered by the use of an LSH table built over the sketches, with the exception of exact mode (--countdict or --set), which use an LSH index built over their bottom-k hashes.\n"\
        "For details on LSH table parameters, see `LSH Options` below.\n"\
        "Top-K (K-Nearest-Neighbor) mode -- \n"\
        "--topk <arg>\tMaximum number of nearest neighbors to list. If <arg> is greater than N - 1, pairwise distances are instead emitted.\n"\
        "Thresholded Mode -- \n"\
        "--similarity-threshold <arg>\tMinimum fraction similarity for inclusion.\n\tIf this is enabled, only pairwise similarities over <arg> will be emitted.\n"\
        "Distance Output Options--\n"\
        "--binary-output\tEmit binary output rather than human-readable.\n\n"\
        "\t For symmetric pairwise, this emits condensed distance matrix in f32\n"\
        "\t For -Q/--qfile usage, this emits a full rectangular distance matrix in of shape (|F|, |Q|)\n"\
        "\t For asymmetric pairwise, this emits a full distance matrix in f32 in row-major storage.\n"\
        "\t For asymmetric pairwise, this emits a flat distance matrix in f32\n"\
        "\t For top-k filtered, this emits a matrix of min(k, |N|) x |N| of IDs and distances\n"\
        "In `dashing2 sketch`, distances are not automatically computed; If set, they are written to <arg>\n"\
        "In `dashing2 cmp`, this defaults to stdout.\n"\
        "--cmpout/--distout/--cmp-outfile\tCompute distances and emit them to <arg>.\n"\
        "\t For similarity-thresholded distances, this emits a compressed-sparse row (CSR) formatted matrix with 64-bit indptr, 32-bit indices, and 32-bit floats for distances\n"\
        "We also provide Greedy Clustering for grouping/de-duplication, using the CD High-Identity with Tolerance (CD-HIT) algorithm, also powered by the LSH table.\n"\
        "Greedy HIT Clustering\n"\
        "This is enabled by a single parameter, --greedy <float (0-1]>, which determines the similarity threshold below which a new cluster is formed.\n"\
        "--greedy <float (0-1]> For greedy clustering by a given similarity threshold; this selects representative sequences or sequence sets.\n"\
        "As this number approaches 1, the number and uniformity of clusters grows.\n"\
        "For human-readable, this emits one line per cluster listing its constituents, ordered by similarity\n"\
        "For machine-readable, this file consists of 2 64-bit integers (nclusters, nsets), followed by (nclusters + 1) 64-bit integers, followed by nsets 64-bit integers, identifying which sets belonged to which clusters.\n"\
        "This is a vector in Compressed-Sparse notation.\n"\
        "Example Python code:\n"\
        "def parsef(fpath, d64=False):\n"\
        "    import numpy as np;data = np.memmap(fpath)\n"\
        "    nclusters, nsets = data[:16].view(np.uint64)[:2]\n"\
        "    endp = 16 + 8 * nclusters\n"\
        "    indptr = data[16:endp].view(np.uint64)\n"\
        "    #  dashing2 uses 32-bit ids\n"\
        "    #  dashing2-64 uses 64-bit ids\n"\
        "    indicesdtype = np.uint64 if d64 else np.uint32\n"\
        "    indices = data[endp:].view(indicesdtype)\n"\
        "    return [indices[start:end] for start, end in zip(indptr[:-1],indptr[1:])]\n"\
        "Runtime Options --\n"\
        "By default, we compare items with full hash function registers; to trade accuracy for speed, these sketches can be compressed before comparisons.\n"\
        "--fastcmp <arg>\tEnable faster comparisons using n-byte signatures rather than full registers. By default, these are logarithmically-compressed\n"\
        "For example, --fastcmp 1 uses byte-sized sketches, with a and b parameters inferred by the data to minimize information loss\n"\
        "<arg> may be 8 (64 bits), 4 (32 bits), 2 (16 bits), 1 (8 bits), or .5 (4 bits)\n"\
        "Results may even be somewhat stabilized by the use of smaller registers.\n"\
        "      Dependent Option:\n"\
        "          --bbit-sigs: truncate to bottom-<arg> bytes of signatures instead of logarithmically-compressed.\n"\
        "For minhash sketches (setsketch, probminhash, and bagminhash), this uses full registers instead of compressed.\n"\
        "--refine-exact\tThis causes the candidate KNN graph to be refined to a final KNN graph using full distances.\tIf using sketches, then full hash registers are used.\nOtherwise, exact k-mer comparison functions are used.\n"\
        "\n\nDistance specifications\n"\
        "The default value emitted is similarity. For MinHash/HLL/SetSketch sketches, this is the fraction of shared registers.\n"\
        "This can be changed to a distance (--mash-distance) for k-mer similarity, where it can be used for hierarchical clustering.\n"\
        "--mash-distance/--poisson-distance\t Emit distances, as estimated by the Poisson model for k-mer distances.\n"\
        "--symmetric-containment\t Use symmetric containment as the distance. e.g., (|A & B| / min(|A|, |B|))\n"\
        "--containment\t Use containment as the distance. e.g., (|A & B| / |A|). This is asymmetric, so you must consider that when deciding the output shape.\n"\
        "--compute-edit-distance\t For edit distance, perform actual edit distance calculations rather than returning the distance in LSH space.\n"\
        "                       \t This means that the LSH index eliminates the quadratic barrier in candidate generation, but they are refined using actual edit distance.\n"\
        "--batch-size [16] \tFor rectangular distance calculation (symmetric all-pairs, asymmetric all-pairs, and query-reference), this batches computation so that memory requirements overlap for better cache-efficiency.\n"\
        "                  \tBy default, this is 16. Increasing the batch size may yield substantial performance improvements, especially is the sketches are rather small.\n"\
        "LSH Options\n"\
        "There are a variety of heuristics in the LSH tables; however, the most important besides sketch size is the number of hash tables used.\n"\
        "--nLSH <int=2>\t\n"\
        "   Aliases: --nlsh\n"\
        "This sets the number of LSH tables. The first 3 tables use register grouping sizes of powers of 2, and subsequent tables use 2 times the index.\n"\
        "If 2 is used (default), these will be of sizes (1, 2), but 4 yields (1, 2, 4, 6) and 5 yields (1, 2, 4, 6, 8) registers per composite key.\n"\
        "Increase this number to pay more memory/time for higher accuracy.\n"\
        "Decrease this number for higher speed and lower accuracy.\n"\
        "This is ignored for exact sketching (--countdict or --set), where a single permutation is generated and a single hash table is used.\n"\
        "\nMetadata Options\n"\
        "If sketching, you can also choose to save k-mers (the IDs corresponding to the k-mer selected), or\n"\
        " and optionally save the counts for these k-mers\n"\
        "This could be used to build inverted indexes (using samples to estimate containment), or for frequency estimation\n"\
        "-s/--save-kmers: Save k-mers. This puts the k-mers saved into .kmer files to correspond with the minhash samples. \n"\
        "If an output path is specified for dashing2 and --save-kmers is enabled, stacked k-mers will be written to <arg>.kmer64, and names will be written to <arg>.kmer.names.txt\n"\
        "This has a 16-byte header containing a 32-bit integer describing the alphabet used, 32 bits describing sketch size, one 32-bit integer for k, and one 32-bit integer for window-length.\n"\
        "This database can be used for dashing2 contain.\n"\
        "-N/--save-kmercounts: Save k-mer counts for sketches. This puts the k-mer counts saved into .kmercounts.f64 files to correspond with the k-mers.\n"\
        "-o/--outfile: sketches are stacked into a single file and written to <arg>\n"\
        "This is the path for the stacked sketches; to set output location, use --cmpout instead. (This is the distance matrix betweek sketches).\n"\


extern size_t MEMSIGTHRESH;

}

#endif /* DASHING2_OPTIONS_H__ */
