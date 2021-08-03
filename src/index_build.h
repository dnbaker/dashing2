#ifndef DASHING2_INDEX_BUILD_H__
#define DASHING2_INDEX_BUILD_H__
#include "src/ssi.h"
#include "src/cmp_main.h"
namespace dashing2 {

using PairT = std::pair<LSHDistType, LSHIDType>;


struct pqueue: public std::priority_queue<PairT> {
    using base_type = std::vector<PairT>;
    base_type &getc() {return this->c;}
    const base_type &getc() const {return this->c;}
    PairT &operator[](size_t i) {return this->c[i];}
    const PairT &operator[](size_t i) const {return this->c[i];}
    void sort() {std::sort(this->c.begin(), this->c.end());}
    PairT &front() {return this->c[0];}
    const PairT &front() const {return this->c[0];}
    PairT &back() {return this->c[this->c.size() - 1];}
    const PairT &back() const {return this->c[this->c.size() - 1];}
    void reserve(size_t n) {this->c.reserve(n);}
    void resize(size_t n) {this->c.resize(n);}
    auto begin() {return this->c.begin();}
    auto begin() const {return this->c.begin();}
    auto end() {return this->c.end();}
    auto end() const {return this->c.end();}
    void erase(typename std::priority_queue<PairT>::container_type::iterator it, typename std::priority_queue<PairT>::container_type::iterator oit);
};


std::vector<pqueue> build_index(SetSketchIndex<LSHIDType, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result, const bool index_compressed=false);

} // namespace dashing2

#endif
