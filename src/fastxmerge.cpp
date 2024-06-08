#include "fastxmerge.h" //added this line, along creating the header file
#include "fastxsketch.h"

namespace dashing2 {
//This merge function only merges multiple sketches into one vector of sketches, it doesn't actually merge sketches. I think this function is needed after the parallel sketching part
SketchingResult SketchingResult::merge(SketchingResult *start, size_t n, const std::vector<std::string> &names=std::vector<std::string>()) {
    DBG_ONLY(std::fprintf(stderr, "About to merge from %p of size %zu, names has size %zu\n", (void *)start, n, names.size());)
    SketchingResult ret;
    if(n == 0) return ret; //no Sketches: return result in SketchingResult directly
    else if(n == 1) { //One sketch: move it to the result and prepend the first name from the names vector to each name in the sketch, then return.
        ret = std::move(*start);
        std::transform(ret.names_.begin(), ret.names_.end(), ret.names_.begin(), [&names](const auto &x) {return names.front() + ":" + x;});
        return ret;
    }

    //Prepare for Merging: Loop through files and collect information
    //ret.nperfile_.resize(total_seq);
    for(size_t i = 0; i < n; ++i) {
        ret.nperfile_.insert(ret.nperfile_.end(), start[i].nperfile_.begin(), start[i].nperfile_.end());
    }
    size_t total_seqs = 0, total_sig_size = 0;
    std::vector<size_t> offsets(n + 1);
    std::vector<size_t> sig_offsets(n + 1);
    for(size_t i = 0; i < n; ++i) {
        const size_t nseqsi = start[i].names_.size();
        const size_t nregs = start[i].signatures_.size();
        assert(start[i].names_.size() == start[i].cardinalities_.size());
        total_seqs += nseqsi;
        total_sig_size += nregs;
        offsets[i + 1] = total_seqs;
        sig_offsets[i + 1] = total_sig_size;
    }
    //Resize the result containers
    ret.names_.resize(total_seqs);
    if(std::any_of(start, start + n, [](auto &x) {return x.sequences_.size();})) {
        seq_resize(ret.sequences_, total_seqs);
    }
    const size_t sketchsz = start->signatures_.size() / start->names_.size();
    if(total_sig_size > 0) {
        ret.signatures_.resize(total_sig_size);
    }
    if(start->kmers_.size()) {
        ret.kmers_.resize(total_seqs * sketchsz);
    }
    ret.cardinalities_.resize(total_seqs);
    if(start->kmercounts_.size()) {
        ret.kmercounts_.resize(total_sig_size);
    }
    const bool seqsz = total_seqs,
               kmercountsz = !start->kmercounts_.empty();
    //Copy data from each sketch -> does the concatenation
    for(size_t i = 0; i < n; ++i) {
        auto &src = start[i];
        assert(src.names_.size() == offsets[i + 1] - offsets[i]);
        const auto ofs = offsets[i], sofs = sig_offsets[i];
        std::string fname;
        if(names.size() > i) fname = names[i].substr(0, names[i].find_first_of(' '));
        // Append filename to sequence names to ensure seq names
        std::transform(src.names_.begin(), src.names_.end(), &ret.names_[ofs], [&fname](const auto &x) {
            return x + ':' + fname;
        });
        std::copy(src.cardinalities_.begin(), src.cardinalities_.end(), &ret.cardinalities_.at(ofs));
        if(seqsz) {
            ret.sequences_.add_set(src.sequences_.begin(), src.sequences_.end());
        }
        if(!start[i].signatures_.empty())
            std::copy(src.signatures_.begin(), src.signatures_.end(), &ret.signatures_.at(sofs));
        if(!start[i].kmers_.empty())
            std::copy(src.kmers_.begin(), src.kmers_.end(), &ret.kmers_[sofs]);
        if(kmercountsz)
            std::copy(src.kmercounts_.begin(), src.kmercounts_.end(), &ret.kmercounts_[sofs]);
    }
    return ret;
}

//creates the filename for the sketchfile, but doesn't create the file yet
//opts: holds infos about sketch, path: base path for filename, iskmer: say wether file is related to k-mers
std::string makedest(Dashing2Options &opts, const std::string &path, bool iskmer) {
    std::string ret(path);
    ret = ret.substr(0, ret.find_first_of(' '));
    if(opts.trim_folder_paths() || opts.outprefix_.size()) {
        ret = trim_folder(path);
        if(opts.outprefix_.size())
            ret = opts.outprefix_ + '/' + ret;
    }
    if(opts.seedseed_ != 0)
        ret += ".seed" + std::to_string(opts.seedseed_);
    if(opts.canonicalize())
        ret += ".rc_canon";
    if(!opts.sp_.unspaced()) {
        ret += opts.sp_.to_string();
    }
    if(opts.kmer_result_ <= FULL_SETSKETCH)
        ret = ret + std::string(".sketchsize") + std::to_string(opts.sketchsize_);
    ret = ret + std::string(".k") + std::to_string(opts.k_);
    if(opts.w_ > opts.k_) {
        ret = ret + std::string(".w") + std::to_string(opts.w_);
    }
    if(opts.count_threshold_ > 0) {
        ret = ret + ".ct_threshold";
        if(std::fmod(opts.count_threshold_, 1.)) ret = ret + std::to_string(opts.count_threshold_);
        else ret = ret + std::to_string(int(opts.count_threshold_));
    }
    if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_EDIT_DISTANCE) {
        ret += '.';
        ret += to_string(opts.ct());
        if(opts.ct() != EXACT_COUNTING)
            ret += std::to_string(opts.cssize_);
    }
    if(opts.sspace_ == SPACE_SET && opts.sketch_compressed()) {
        char buf[256];
        auto l = std::sprintf(buf, ".a=%0.16Lg.b=%0.16Lg.fd=%0.16Lg", opts.compressed_a_, opts.compressed_b_, static_cast<long double>(opts.fd_level_));
        ret += std::string(buf, l);
    }
    ret += ".";
    if(opts.kmer_result_ <= FULL_SETSKETCH)
        ret += to_string(opts.sspace_);
    else {
        auto ks = opts.kmer_result_;
        if(iskmer && ks == FULL_MMER_COUNTDICT) {
            ks = FULL_MMER_SET;
        }
        ret += to_string(ks);
    }
    ret = ret + "." + bns::to_string(opts.rht_) + to_suffix(opts);
    DBG_ONLY(std::fprintf(stderr, "Source %s->%s\n", path.data(), ret.data());)
    return ret;
}


}
