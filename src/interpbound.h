#ifndef INTERPBOUND_H__
#define INTERPBOUND_H__

namespace interp {

template<typename IT, typename KeyT>
std::pair<IT, bool> search(const IT beg, const IT end, KeyT key) {
    assert(std::is_sorted(beg, end));
    assert(beg < end);
    if(beg == end) return {end, false};
    using CT = typename std::conditional<(sizeof(KeyT) > 8), long double, double>::type;
    //std::fprintf(stderr, "[%s] span from %zu to %zu\n", __PRETTY_FUNCTION__, *beg, *(end - 1));
    KeyT max = *(end - 1);
    KeyT min = *beg;
    //std::fprintf(stderr, "min: %zu. max: %zu\n", min, max);
    assert(min < max);
    if(key > max) return {end, false};
    if(key < min) return {beg, false};
    IT upper = end - 1;
    IT lower = beg;
    //if(key == *lower) return {lower, true};
    //if(key == *upper) return {upper, true};
    size_t span = upper - lower;
    IT cid = lower + std::ptrdiff_t(CT(key - min) / CT(max - min) * (CT(upper - lower) + 1));
    if(cid > upper) cid = upper;
    //std::fprintf(stderr, "key %zu vs min %zu vs max %zu, yielding %zu/%zu\n", key, min, max, cid - beg, end - beg);
    DBG_ONLY(size_t iternum = 0;)
    for(;;) {
        DBG_ONLY(std::fprintf(stderr, "Iter %zu, %zd/%zu. key: %zu vs current point %zu\n", ++iternum, cid - beg, end - beg, (key), (*cid));)
        assert(cid <= end);
        assert(cid >= beg);
        auto value = *cid;
        if(*cid == key) return {cid, true};
        if(*cid > key) {
            //std::fprintf(stderr, "Lowering upper bound\n");
            upper = cid - 1;
            max = value;
        } else {
            //(*cid < key)
            lower = std::min(cid + 1, end);
            min = value;
        }
        assert(key >= min || key <= max);
        CT rat = CT(key - min) / CT(max - min);
        size_t range = upper - lower;
        std::ptrdiff_t diff = std::ptrdiff_t(rat * (range));
        if(upper < lower) {
            //std::fprintf(stderr, "Upper is lower by %zd...\n", lower - upper);
            return {cid, false};
        }
        cid = lower + diff;
    }
}
template<typename IT, typename KeyT>
IT find(IT beg, IT end, KeyT key) {
    auto [it, v] = search(beg, end, key);
    return v ? it: end;
}

}

#endif
