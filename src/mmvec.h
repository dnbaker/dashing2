#ifndef MMVEC_H__
#define MMVEC_H__
#include <mio.hpp>
#include <cstring>
#include <cassert>
#include <string>
#include <stdexcept>

namespace mm {
using namespace std::literals::string_literals;

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
    std::vector<T> ram_;
public:
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
            std::FILE *ifp;
            if((ifp = std::fopen(path(), "wb")) == nullptr)
                throw std::runtime_error(std::string("Failed to open ") + path_);
            struct stat st;
            if(::fstat(::fileno(ifp), &st)) throw std::runtime_error("Failed to fstat");
            offset_ = off == size_t(-1) ? size_t(st.st_size): off;
            if(std::fclose(ifp)) throw std::runtime_error("Error in closing ifp\n");
        }
        std::error_code ec;
        if(initial_size) {
            capacity_ = size_ = initial_size;
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
    vector(vector &o): offset_(o.offset_), capacity_(o.capacity_), size_(o.size_), path_(o.path_), ms_(std::move(o.ms_)) {
        o.offset_ = o.capacity_ = o.size_ = 0;
        o.path_.clear();
    }
    vector(const vector &o) = delete;
    vector(std::string inpath, size_t off=size_t(-1), size_t initial_size=0, T init=T(0)): capacity_(0), size_(0), path_(inpath) {
        std::FILE *ifp;
        if((ifp = std::fopen(path(), "wb")) == nullptr)
            throw std::runtime_error(std::string("Failed to open ") + path_);
        struct stat st;
        if(::fstat(::fileno(ifp), &st)) throw std::runtime_error("Failed to fstat");
        if(off == size_t(-1))
            off = st.st_size;
        offset_ = off;
        if(std::fclose(ifp)) throw std::runtime_error("Error in closing ifp\n");
        if(initial_size) while(size_ < initial_size) push_back(init);
#ifndef NDEBUG
        std::fprintf(stderr, "Construction finished. Total number of items: %zu\n", size_);
#endif
    }
    template<typename...Args>
    T &emplace_back(Args &&...args) {
        if(size() == 0u || size() == capacity_) {
            size_t newcap = std::max(static_cast<size_t>(1), size() << 1);
            std::fprintf(stderr, "Reached cap %zu, now resizing to %zu\n", capacity_, newcap);
            reserve(newcap);
        }
        std::fprintf(stderr, "Current size: %zu. Capacity: %zu. is kept in RAM: %d\n", size_, capacity_, using_ram());
        T *ptr = data() + size_++;
        new(ptr) T(std::forward<Args>(args)...);
        return *ptr;
    }
    T &push_back(T x) {return emplace_back(std::move(x));}
    ~vector() noexcept(false) {
        // shrink to fit, then flush to disk if necessary
        if(path_.empty()) return;
        if(capacity_ < countthreshold()) {
            ram_.resize(size_);
            ram_.shrink_to_fit();
            capacity_ = size_;
            if(!path_.empty()) {
                size_t diskspace = offset_ + capacity_ * sizeof(T);
                if(::truncate(path(), diskspace))
                    throw std::runtime_error("Failed to resize "s + path_ + "via ::truncate");
#ifndef NDEBUG
                std::fprintf(stderr, "Capacity is < threshold, so we need to flush the data in RAM to the file.\n");
#endif
                std::error_code ec;
                ms_.map(path_, offset_, sizeof(T) * capacity_, ec);
#ifndef NDEBUG
                if(ec) std::fprintf(stderr, "Error mmaping %s\n", ec.message().data());
#endif
                std::copy(ram_.begin(), ram_.end(), (T *)ms_.data());
#ifndef NDEBUG
                std::fprintf(stderr, "Copied from RAM file to disk.\n");
#endif
            }
        } else if(capacity_ > 0 && capacity_ > size_) {
            const bool wasopen = ms_.is_open();
            if(wasopen) ms_.unmap();
#ifndef NDEBUG
            std::fprintf(stderr, "Shrinking %zu to %zu\n", capacity_, size_);
#endif
            if(::truncate(path(), offset_ + (size_ * sizeof(T))))
                throw std::runtime_error("Failed to truncate to shrink_to_fit");
#ifndef NDEBUG
            std::fprintf(stderr, "Size, cap at end: %zu, %zu. Path = %s (%p)\n", size_, capacity_, path(), (void *)this);
#endif
            if(wasopen) {
                std::error_code ec;
                ms_.map(path(), offset_, size_ * sizeof(T), ec);
                if(ec) throw std::runtime_error(ec.message());
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
        if(::truncate(path(), newnbytes))
            throw std::runtime_error("Failed to resize"s + path_ + " via ::truncate");
        std::error_code ec;
        ms_.map(path_, offset_, sizeof(T) * newcap, ec);
        if(ec) throw std::runtime_error(ec.message());
        if(capacity_ < ct && newcap >= ct) {
#ifndef NDEBUG
            std::fprintf(stderr, "Crossing from %zu to %zu passes count threshold %zu, so we are writing the vector from RAM to disk and switching to mmap.\n", capacity_, newcap, ct);
#endif
            std::copy(ram_.begin(), ram_.end(), (T *)ms_.data());
            std::vector<T> tmp(std::move(ram_));
            assert(ram_.empty());
        }
        capacity_ = newcap;
    }
    void shrink_to_fit() {
    }
};

} // namespace mm
#endif /* MMVEC_H__ */
