#include "enums.h"

namespace dashing2 {

std::string to_string(KmerSketchResultType t) {
    if(t == ONE_PERM) {return "OnePermutationSetSketch";}
    if(t == FULL_SETSKETCH) return "FullSetSketch";
    if(t == FULL_MMER_SET) return "FullMmerSet";
    if(t == FULL_MMER_SEQUENCE) return "FullMmerSequence";
    return "FullMmerCountdict";
}
std::string to_string(DataType dt) {
    if(dt == FASTX) return "Fastx";
    if(dt == BED) return "BED";
    if(dt == BIGWIG) return "BigWig";
    if(dt == LEAFCUTTER) return "LeafCutter";
    return "Unknown";
}


}
