#ifndef TEMPSEQS_H__
#define TEMPSEQS_H__
#include <cstdio>
#include <stdexcept>
#include <filesystem>

namespace tmpseq {

struct FileDeleter {
    void operator()(std::FILE* const ptr) const noexcept {
        if(ptr) {
            std::fclose(ptr);
        }
    }
};

class Seqs {
    std::string path_;
    std::vector<int64_t> offsets_;
    std::unique_ptr<std::FILE, FileDeleter> file_;

    static std::string make_path(int64_t seed = 0) {
        std::filesystem::path base = std::filesystem::temp_directory_path();
        std::string randstr;
        for(int i = 0; i < 10; ++i) {
            randstr += (static_cast<char>(seed % 26) + 'a');
            seed += 0x33b74ea7d2a06fd7;
            seed *= 0x7262c70021180663;
            seed ^= (seed >> 16);
        }
        base /= randstr;
        return static_cast<std::string>(base);
    }
    static std::string make_new_file(int64_t seed=0) {
        std::string ret;
        do {
            ret = make_path(seed++);
        } while(std::filesystem::exists(ret));
        return ret;
    }
public:
    Seqs(): path_(make_new_file()), offsets_{0}, file_{std::fopen(path_.data(), "wb")} {}
    Seqs(Seqs&& o) = default;
    Seqs& operator=(Seqs&& o) = default;
    Seqs(Seqs& o) = delete;
    int64_t add_sequence(std::string_view seq) {
        const int64_t ret = offsets_.back();
        const int64_t seq_len = seq.size();
        const int64_t written = std::fwrite(seq.data(), 1, seq_len, file_.get());
        if(written != seq_len) {
            throw std::runtime_error(std::string("ERROR: Failed to append sequence of size ") + std::to_string(seq.size()) + " to temporary file. Out of disk space?");
        }
        offsets_.push_back(ret + seq.size());
        std::fflush(file_.get());
        return ret;
    }
    void emplace_back(const char* ptr, size_t n) {
        add_sequence(std::string_view(ptr, n));
    }
    void push_back(std::string_view seq) {
        add_sequence(seq);
    }
    std::string operator[](const int64_t idx) const {
        thread_local std::unique_ptr<std::FILE, FileDeleter> reader;
        if(!reader) {
            std::FILE *fp = std::fopen(path_.data(), "rb");
            if(fp == nullptr) throw std::runtime_error(std::string("Failed to open file at ") + path_ + " for reading.");
            reader.reset(fp);
        }
        const int64_t offset = offsets_[idx];
        const int64_t end = offsets_[idx + 1];
        std::fseek(reader.get(), offset, SEEK_SET);
        const int64_t ret_size = end - offset;
        std::string ret(ret_size, '\0');
        const int64_t num_read = std::fread(ret.data(), 1, ret.size(), reader.get());
        if(num_read != ret_size) {
            throw std::runtime_error(std::string("ERROR: Failed to read sequence ") + std::to_string(idx) + '/' + std::to_string(size()) + " total. Read " + std::to_string(num_read) + " and expected " + std::to_string(ret.size()));
        }
        return ret;
    }
    int64_t size() const noexcept {
        return (offsets_.size()) - 1;
    }
};

}

#endif
