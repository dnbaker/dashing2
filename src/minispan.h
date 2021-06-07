#ifndef DASHING2_MINISPAN_H__
#define DASHING2_MINISPAN_H__
#include <cstdlib>

namespace dashing2 {
using std::size_t;

template<typename T>
struct minispan {
    T *ptr_;
    size_t n_;
    minispan(T *ptr, size_t n): ptr_(ptr), n_(n) {}
    minispan(const T *ptr, size_t n): ptr_(const_cast<T *>(ptr)), n_(n) {}
    size_t size() const {return n_;}
    T *begin() {return ptr_;}
    const T *begin() const {return ptr_;}
    T *end() {return ptr_ + n_;}
    const T *end() const {return ptr_ + n_;}
    T *data() {return ptr_;}
    const T *data() const {return ptr_;}
    T &operator[](size_t idx) {return ptr_[idx];}
    const T &operator[](size_t idx) const {return ptr_[idx];}
};

}

#endif /* DASHING2_MINISPAN_H__ */
