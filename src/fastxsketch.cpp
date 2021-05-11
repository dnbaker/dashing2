#include "d2.h"
#include <optional>
namespace dashing2 {

struct FastxSketchingResult {
    std::vector<std::string> names_;
    std::vector<uint32_t> nperfile_; // This is either empty (in which case each filename/row has its own sketch)
                                     // Or, this contains a list indicating the number of sketches created for each file/line
    std::vector<RegT> signatures_;
};

struct SequenceSketch {
    std::optional<Counter> ctr; // If counts were generated, rather than
    std::optional<std::vector<RegT>> sigs;
};

FastxSketchingResult fastx2sketch(const ParseOptions &opts, std::vector<std::string> &paths) {
    FastxSketchingResult ret;
    return ret;
}

} // dashing2
