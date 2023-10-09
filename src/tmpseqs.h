#ifndef TEMPSEQS_H__
#define TEMPSEQS_H__
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <string>
#include <variant>
#include <vector>
#include <mio.hpp>

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
    Seqs(bool in_memory=false): path_(make_new_file(in_memory)), offsets_{0}, file_{std::fopen(path_.data(), "wb")} {}
    Seqs(Seqs&& o) = default;
    Seqs& operator=(Seqs&& o) = default;
    Seqs(Seqs& o) = delete;
    ~Seqs() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
        if(ec) {
            std::fprintf(stderr, "Failed to delete temporary file %s. You may need to delete this manually.\n", path_.data());
        }
    }
    void swap_to_ram() const noexcept {
    }
    void free_if_possible(const int64_t _=0) const noexcept {
        if(_) {
            std::fprintf(stderr, "Warning: cannot free a ram-backed database.\n");
        }
    }
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
    void emplace_back(std::string_view seq) {
        add_sequence(seq);
    }
    template<typename It>
    void add_set(It beg, const It end) {
        while(beg != end) {
            emplace_back(*beg);
            ++beg;
        }
    }
#if 1
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
#else
    std::string operator[](const int64_t idx) const {
        static std::unique_ptr<mio::mmap_sink> reader;
        static std::mutex mutex;
        if(!reader) {
            std::lock_guard<std::mutex> lock(mutex);
            if(!reader) {
                reader = std::make_unique<mio::mmap_sink>(path_);
            }
        }
        return std::string(reader->data() + offsets_[idx], offsets_[idx + 1] - offsets_[idx]);
    }
#endif
    int64_t size() const noexcept {
        return (offsets_.size()) - 1;
    }
    struct Iterator {
        const Seqs& container;
        int64_t idx;
        std::string_view operator*() const {
            return container[idx];
        }
        Iterator& operator++() {
            ++idx;
            return *this;
        }
        Iterator operator++(int) {
            Iterator ret{*this};
            ++idx;
            return ret;
        }
        bool operator==(const Iterator& other) {
            return idx == other.idx;
        }
        bool operator!=(const Iterator& other) {
            return idx != other.idx;
        }
    };
    Iterator begin() {
        return {*this, static_cast<int64_t>(0)};
    }
    Iterator end() {
        return {*this, size()};
    }
};

struct MemoryOrRAMSequences {
    using V = std::variant<std::vector<std::string>, Seqs>;
    V core_;
    MemoryOrRAMSequences(const bool inRam): core_{inRam ? V(std::vector<std::string>{}): V(Seqs{})} {}
    struct Iterator {
        const MemoryOrRAMSequences& container;
        int64_t idx;
        std::string operator*() const {
            return container[idx];
        }
        Iterator& operator++() {
            ++idx;
            return *this;
        }
        Iterator operator++(int) {
            Iterator ret{*this};
            ++idx;
            return ret;
        }
        bool operator==(const Iterator& other) {
            return idx == other.idx;
        }
        bool operator!=(const Iterator& other) {
            return idx != other.idx;
        }
    };
    std::string operator[](const int64_t idx) const {
        return std::visit([idx](const auto& x) -> std::string {
            return std::string(x.operator[](idx));
        }, core_);
    }
    Iterator begin() {
        return {*this, static_cast<int64_t>(0)};
    }
    Iterator end() {
        return {*this, size()};
    }
    void emplace_back(const char* seq, const size_t len) {
        emplace_back(std::string_view(seq, len));
    }
    void emplace_back(std::string_view seq) {
        std::visit([seq](auto& x) -> void {
            x.emplace_back(seq);
        }, core_);
    }
    int64_t size() const noexcept {
        return std::visit([](const auto& x) -> int64_t {return x.size();}, core_);
    }
    template<typename It>
    void add_set(It beg, const It end) {
        while(beg != end) {
            emplace_back(*beg);
            ++beg;
        }
    }
    void free_if_possible(const int64_t idx) {
        if(core_.index() == 0) {
            std::string tmp;
            std::swap(tmp, std::get<0>(core_)[idx]);
        }
    }
    void free_if_possible() {
        if(core_.index() == 0) {
            std::vector<std::string> tmp;
            std::swap(tmp, std::get<0>(core_));
        }
    }
    void swap_to_ram() {
        if(size() > 0) {
            throw std::runtime_error("Cannot swap variant type after adding an item.");
        }
        if(core_.index() == 1) {
            core_ = V(std::vector<std::string>{});
        }
    }
};

} // tmpseq

#endif
