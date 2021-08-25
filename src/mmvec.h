#ifndef MMVEC_H__
#define MMVEC_H__
#include <mio.hpp>
#include <cstring>
#include <string>
#include <stdexcept>

namespace mm {
using namespace std::literals::string_literals;

template<typename T>
class vector {
    // WARNING: this doesn't handle constructors or destructors. This is only for POD data!
    size_t offset;
    size_t capacity_;
    size_t size_;
    std::string path_;
    mio::mmap_sink ms_;
public:
    size_t capacity() const {return capacity_;}
    size_t size() const {return size_;}
    bool empty() const {return size_ == 0u;}
    int fd() const {return ms_.file_handle();}
    const void *raw_data() const {return static_cast<const void *>(ms_.data());}
    void *raw_data() {return static_cast<void *>(ms_.data());}
    const T *data() const {
        return static_cast<const T *>(raw_data());
    }
    T *data() {
        return static_cast<T *>(raw_data());
    }
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
    void assign(std::string path, size_t off=size_t(-1), size_t initial_size=0, T init=T(0)) {
        if(ms_.is_open()) ms_.unmap();
        path_ = path;
        std::FILE *ifp;
        if((ifp = std::fopen(path_.data(), "wb")) == nullptr)
            throw std::runtime_error(std::string("Failed to open ") + path_);
        struct stat st;
        if(::fstat(::fileno(ifp), &st)) throw std::runtime_error("Failed to fstat");
        if(off == size_t(-1))
            off = st.st_size;
        offset = off;
        if(std::fclose(ifp)) throw std::runtime_error("Error in closing ifp\n");
        std::error_code ec;
        if(initial_size) {
            capacity_ = size_ = initial_size;
            if(::truncate(path_.data(), sizeof(T) * initial_size + offset)) throw std::runtime_error("Failed to resize file.");
            ms_.map(path_, offset, sizeof(T) * initial_size, ec);
            if(ec) throw std::runtime_error(ec.message() + "; failed to mmap");
            for(size_t i = 0; i < initial_size; new(data() + i++) T(init));
        } else {
            ::truncate(path_.data(), offset);
        }
#ifndef NDEBUG
        std::fprintf(stderr, "Construction finished. Total number of items: %zu\n", size_);
#endif
    }
    vector(): offset(-1), capacity_(0), size_(0) {
    }
    vector(vector &&o): vector((vector &)o) {}
    vector(vector &o): offset(o.offset), capacity_(o.capacity_), size_(o.size_), path_(o.path_), ms_(std::move(o.ms_)) {
        o.offset = o.capacity_ = o.size_ = 0;
        o.path_.clear();
    }
    vector(const vector &o) = delete;
    vector(std::string path, size_t off=size_t(-1), size_t initial_size=0, T init=T(0)): size_(0), capacity_(0), path_(path) {
        std::FILE *ifp;
        if((ifp = std::fopen(path_.data(), "wb")) == nullptr)
            throw std::runtime_error(std::string("Failed to open ") + path_);
        struct stat st;
        if(::fstat(::fileno(ifp), &st)) throw std::runtime_error("Failed to fstat");
        if(off == size_t(-1))
            off = st.st_size;
        offset = off;
        if(std::fclose(ifp)) throw std::runtime_error("Error in closing ifp\n");
        if(initial_size) {
            std::error_code ec;
            capacity_ = size_ = initial_size;
            if(::truncate(path_.data(), sizeof(T) * initial_size + offset)) throw std::runtime_error("Failed to resize file.");
            ms_.map(path_, offset, sizeof(T) * initial_size, ec);
            if(ec) throw std::runtime_error(ec.message() + "; failed to mmap");
            for(size_t i = 0; i < initial_size; new(data() + i++) T(init));
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
    T &push_back(T &&x) {
        return emplace_back(x);
    }
    ~vector() {
        if(ms_.is_open()) {
            ms_.unmap();
        }
        shrink_to_fit();
        offset = capacity_ = size_ = 0;
        path_ = "";
    }
    void resize(size_t newsz, T init=T(0)) {
        if(newsz > size_) {
            reserve(newsz);
            std::fill_n(data() + size_, newsz - size_, init);
        }
        size_ = newsz;
    }
    void reserve(size_t newsz) {
        if(newsz <= capacity_) return;
        if(ms_.is_open()) ms_.unmap();
        size_t newnbytes = offset + newsz * sizeof(T);
        if(::truncate(path_.data(), newnbytes))
            throw std::runtime_error("Failed to resize"s + path_ + " via ::truncate");
        std::error_code ec;
        ms_.map(path_, offset, sizeof(T) * newsz, ec);
        if(ec) throw std::runtime_error(ec.message());
        capacity_ = newsz;
    }
    void shrink_to_fit() {
        if(capacity_ > 0 && capacity_ > size_) {
#ifndef NDEBUG
            std::fprintf(stderr, "Shrinking %zu to %zu\n", capacity_, size_);
#endif
            if(::truncate(path_.data(), offset + (size_ * sizeof(T)))) throw std::runtime_error("Failed to truncate to shrink_to_fit");
#ifndef NDEBUG
            std::fprintf(stderr, "Size, cap at end: %zu, %zu. Path = %s (%p)\n", size_, capacity_, path_.data(), (void *)this);
#endif
            capacity_ = size_;
        }
    }
};

} // namespace mm
#endif /* MMVEC_H__ */
