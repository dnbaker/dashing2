#ifndef MMVEC_H__
#define MMVEC_H__
#include <mio.hpp>
#include <cstring>
#include <cassert>
#include <string>
#include <stdexcept>
#include "sketch/sseutil.h"

namespace mm {
using namespace std::literals::string_literals;

class AnonMMapper {
    void *data_;
    size_t n_;
public:
    void *data() {return data_;}
    const void *data() const {return data_;}
    static void *allocate(size_t nbytes) {
        return ::mmap(nullptr, nbytes, PROT_READ | PROT_WRITE, MAP_SHARED | MMAP_HUGE_FLAGS | MAP_ANON, -1, 0);
    }
    AnonMMapper(size_t nbytes): n_(nbytes) {
        if(nbytes == 0) {
            return;
        }
        void *tmp = allocate(nbytes);
        if(tmp == reinterpret_cast<void *>(0xffffffffffffffff)) {
            perror("Failed to perform mmap call; throwing std::bad_alloc");
            throw std::bad_alloc();
        }
        data_ = tmp;
    }
    ~AnonMMapper() {
        if(data_) {
            ::munmap(data_, n_);
        }
    }
    AnonMMapper(AnonMMapper &&o): data_(o.data_), n_(o.n_) {
        o.data_ = nullptr;
        o.n_ = 0u;
    }
};

template<typename T>
class vector {
    static_assert(std::is_trivial_v<T>, "T must be trivial");
    static_assert(std::is_standard_layout_v<T>, "T must be standard_layout");
    // WARNING: this doesn't handle constructors or destructors. This is only for POD data!
    size_t offset_;
    size_t capacity_;
    size_t size_;
    std::string path_;
    mio::mmap_sink ms_;
    size_t memthreshold_ = 20ull << 30;
    using VecT = std::vector<T, sse::AlignedAllocator<T>>;
    VecT ram_;
    // TODO: add in additional option for backing to use anonymous mmap'd data
    // using fd = -1, which makes a temporary file which is cleared when no longer needed.
public:
    using value_type = T;
    using reference_type = T &;
    using pointer_type = T *;
    using const_reference_type = T &;
    using const_pointer_type = T *;
    using different_type = std::ptrdiff_t;
    const char *path() const {return path_.data();}
    size_t capacity() const {return capacity_;}
    size_t size() const {return size_;}
    bool empty() const {return size_ == 0u;}
    bool using_ram() const {
        return path_.empty() || capacity_ < countthreshold();
    }
    int fd() const {return ms_.file_handle();}
    const void *raw_data() const {
        assert(ms_.is_open());
        return static_cast<const void *>(ms_.data());
    }
    T &back() {return data()[size_ - 1];}
    const T &back() const {return data()[size_ - 1];}
    void *raw_data() {
        assert(ms_.is_open());
        return static_cast<void *>(ms_.data());
    }
    const T *data() const {
        if(using_ram())
            return ram_.data();
        assert(ms_.is_open());
        return static_cast<const T *>(raw_data());
    }
    T *data() {
        if(using_ram())
            return ram_.data();
        assert(ms_.is_open());
        return static_cast<T *>(raw_data());
    }
    size_t countthreshold() const {
        return memthreshold_ / sizeof(T);
    }
    size_t countthreshold(size_t nitems) {
        memthreshold(nitems * sizeof(T));
        return nitems;
    }
    size_t memthreshold(size_t nb, char type=0) {
        switch(type | 32) {
            case 't': nb <<= 40; break;
            case 'g': nb <<= 30; break;
            case 'm': nb <<= 20; break;
            case 'k': nb <<= 10; break;
        }
        memthreshold_ = nb;
        return memthreshold_;
    }
    size_t memthreshold() const {return memthreshold_;}
    size_t offset() const {return offset_ != size_t(-1) ? offset_: size_t(0);}
    T *begin() {return data();}
    T *end() {return data() + size_;}
    T &operator[](size_t i) {return data()[i];}
    const T &operator[](size_t i) const {return data()[i];}
    T &at(size_t i) {
        if(__builtin_expect(i >= size_, 0)) throw std::out_of_range(std::to_string(i) + " out of range for " + std::to_string(size_));
        return data()[i];
    }
    const T &at(size_t i) const {
        if(__builtin_expect(i >= size_, 0)) throw std::out_of_range(std::to_string(i) + " out of range for " + std::to_string(size_));
        return data()[i];
    }
    vector &operator=(const vector &o) = delete;
    vector &operator=(vector &&o) {
        static constexpr size_t n = sizeof(*this);
        std::swap_ranges((uint8_t *)this, (uint8_t *)this + n, (uint8_t *)&o);
        return *this;
    }
    void assign(std::string inpath="", size_t off=size_t(-1), size_t initial_size=0, T init=T(0)) {
        if(ms_.is_open()) ms_.unmap();
        path_ = inpath;
        if(path_.size()) {
            struct stat st;
            if(::stat(path(), &st)) {
                std::FILE *fp = std::fopen(path(), "wb");
                if(fp == 0) {
                    throw std::runtime_error("Failed to open path "s + path());
                }
                const int fd = ::fileno(fp);
                if(::fstat(fd, &st) < 0) throw std::runtime_error("Failed to stat "s + path());
                std::fclose(fp);
            }
            if(off == size_t(-1)) {
                offset_ = static_cast<size_t>(st.st_size);
                capacity_ = size_ = 0;
            } else {
                offset_ = off;
                if(off >= static_cast<size_t>(st.st_size)) {
                    capacity_ = size_ = 0;
                } else {
                    assert((static_cast<size_t>(st.st_size) - off) % sizeof(T) == 0u);
                    size_ = capacity_ = (static_cast<size_t>(st.st_size) - off) / sizeof(T);
                }
                if(::truncate(path(), size_ * sizeof(T) + offset_)) {
                    perror("Failed to truncate");
                    throw std::runtime_error("Failed to truncate from size "s + std::to_string(st.st_size) + " to size " + std::to_string(size_ * sizeof(T) + offset_));
                }
                countthreshold(size_);
            }
            std::error_code ec;
            ms_.map(path_, offset_, sizeof(T) * capacity_, ec);
            if(ec) {
                perror("Failure");
                throw std::runtime_error("Failed to map path "s + path_);
            }
        } else {
            capacity_ = size_ = initial_size;
        }
        if(initial_size) {
            reserve(initial_size);
            while(size_ < initial_size) push_back(init);
        } else {
            ::truncate(path(), offset_);
        }
#ifndef NDEBUG
        std::fprintf(stderr, "Construction finished. Total number of items: %zu\n", size_);
#endif
    }
    vector(): offset_(-1), capacity_(0), size_(0) {
    }
    vector(vector &&o): vector((vector &)o) {}
    vector(vector &o): offset_(o.offset_), capacity_(o.capacity_), size_(o.size_), path_(o.path_), ms_(std::move(o.ms_)), memthreshold_(o.memthreshold_), ram_(std::move(ram_)) {
        o.offset_ = o.capacity_ = o.size_ = 0;
        o.path_.clear();
    }
    vector(const vector &o) = delete;
    vector(std::string inpath, size_t off=size_t(-1), size_t initial_size=0, T init=T(0)): capacity_(0), size_(0), path_(inpath) {
        struct stat st;
        if(inpath.size()) {
            const char *const ipp = inpath.data();
            if(::stat(ipp, &st)) {
                std::FILE *ifp = std::fopen(ipp, "wb");
                if(0 == ifp || ::fstat(::fileno(ifp), &st))
                    throw std::runtime_error("Failed to create file at "s +path);
                std::fclose(ifp);
            }
            offset_ = off != size_t(-1) ? off: size_t(st.st_size);
        }
        if(initial_size) while(size_ < initial_size) push_back(init);
        else {
            size_ = (static_cast<size_t>(st.st_size) - offset_) / sizeof(T);
            assert((static_cast<size_t>(st.st_size) - offset_ ) % sizeof(T) == 0u);
        }
#ifndef NDEBUG
        std::fprintf(stderr, "Construction finished. Total number of items: %zu\n", size_);
#endif
    }
    template<typename...Args>
    T &emplace_back(Args &&...args) {
        if(size() == 0u || size() == capacity_) {
            size_t newcap = std::max(static_cast<size_t>(1), size() << 1);
            reserve(newcap);
        }
        T *ptr = data() + size_++;
        new(ptr) T(std::forward<Args>(args)...);
        return *ptr;
    }
    T &push_back(T x) {return emplace_back(std::move(x));}
    static size_t fsize(const std::string &path) {
        struct stat st;
        if(::stat(path.data(), &st)) throw std::runtime_error("Failed to stat "s + path + " for fsize");
        return st.st_size;
    }
    static size_t fsize(std::FILE *fp) {
        struct stat st;
        if(::fstat(::fileno(fp), &st)) throw std::runtime_error("Failed to fstat for ptr"s + std::to_string(reinterpret_cast<uint64_t>(fp)));
        return st.st_size;
    }
    static size_t fsize(int fd) {
        struct stat st;
        if(::fstat(fd, &st)) throw std::runtime_error("Failed to fsize");
        return st.st_size;
    }
    ~vector() {
        // shrink to fit, then flush to disk if necessary
        if(path_.empty()) return;
        if(ram_.size() && capacity_ < countthreshold()) {
            ram_.resize(size_);
            ram_.shrink_to_fit();
            capacity_ = size_;
            if(::truncate(path(), offset_ + capacity_ * sizeof(T))) {
                struct stat st;
                ::stat(path(), &st);
                const std::string msg = ("Failed to resize "s + path_ + "via ::truncate to size " + std::to_string(offset_ + capacity_ * sizeof(T)) + " from " + std::to_string(st.st_size) + __FILE__ + ", " + __PRETTY_FUNCTION__);
                std::fprintf(stderr, "%s\n", msg.data());
                std::exit(EXIT_FAILURE);
            }
#ifndef NDEBUG
            std::fprintf(stderr, "Capacity is < threshold, so we need to flush the data in RAM to the file.\n");
#endif
            std::error_code ec;
            ms_.map(path_, offset_, sizeof(T) * capacity_, ec);
            if(ec) {
                std::fprintf(stderr, "Failed to map resized file before flushing to disk. System error: %s\n", std::system_error(ec).what());
                std::exit(EXIT_FAILURE);
            }
            //std::fprintf(stderr, "Flushing to file with %zu as offset.\n", offset_);
            std::copy(ram_.begin(), ram_.end(), (T *)ms_.data());
        } else if(capacity_ > 0 && capacity_ > size_) {
            const bool wasopen = ms_.is_open();
            if(wasopen) ms_.unmap();
            if(::truncate(path(), offset_ + (size_ * sizeof(T)))) {
                std::fprintf(stderr, "Failed to truncate for shrink_to_fit in mmvec destructor.\n");
                std::exit(EXIT_FAILURE);
            }
            if(wasopen) {
                std::error_code ec;
                ms_.map(path(), offset_, size_ * sizeof(T), ec);
                if(ec) {
                    std::fprintf(stderr, "Failed to map resized file before flushing to disk. System error: %s\n", std::system_error(ec).what());
                    std::exit(EXIT_FAILURE);
                }
            }
            capacity_ = size_;
        }
        // Then clean up
        if(ms_.is_open()) ms_.unmap();
        offset_ = capacity_ = size_ = 0;
        path_ = "";
    }
    void resize(size_t newsz, T init=T(0)) {
        if(newsz > size_) {
            reserve(newsz);
            std::fill_n(data() + size_, newsz - size_, init);
        } // else if(newsz < size_) std::destroy_n(data() + newsz, size_ - newsz); (this would be needed for non-trivial types.)
        size_ = newsz;
    }
    void reserve(size_t newcap) {
        if(newcap <= capacity_) return;
        const size_t ct = countthreshold();
        if(path_.empty() || newcap < ct) {
            ram_.resize(newcap);
            capacity_ = ram_.size();
            return;
        }
        if(ms_.is_open()) ms_.unmap();
        size_t newnbytes = offset_ + newcap * sizeof(T);
        if(path_.empty())
            throw std::runtime_error("Failed to spill to disk, as path is empty");
        if(::truncate(path(), newnbytes)) {
            struct stat st;
            ::stat(path(), &st);
            throw std::runtime_error("Failed to resize"s + path_ + " via ::truncate");
        }
        std::error_code ec;
        ms_.map(path_, offset_, sizeof(T) * newcap, ec);
        if(ec) {
            perror("Don't have permission to map.\n");
            throw std::runtime_error(ec.message());
        }
        if(capacity_ <= ct && newcap >= ct && ram_.size()) {
#ifndef NDEBUG
            std::fprintf(stderr, "Crossing from %zu to %zu passes count threshold %zu, so we are writing the vector from RAM to disk and switching to mmap.\n", capacity_, newcap, ct);
#endif
            std::copy(ram_.begin(), ram_.end(), (T *)ms_.data());
            VecT tmp(std::move(ram_));
            assert(ram_.empty());
        }
        capacity_ = newcap;
    }
    void clear() {
        if(using_ram()) {
            ram_.clear();
        }
        size_ = 0;
    }
};

} // namespace mm
#endif /* MMVEC_H__ */
