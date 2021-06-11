#include "src/oph.h"
#include "sketch/setsketch.h"
using namespace dashing2;

int main() {
    for(const int isz: {128, 1024, 8192, 16384}) {
        for(const size_t exp: {100000, 10000000, 500000000}) {
#if 1
    LazyOnePermSetSketch<uint64_t> l(isz);
#else
    sketch::setsketch::OPCSetSketch<double> l(isz);
#endif
    for(size_t i = 0; i < exp; ++i) l.update(i);
    auto p = l.data();
    auto c = l.getcard();
    double ocard2 = isz / std::accumulate(p, p + isz, 0.);
    double err = std::abs(c - exp);
    double ferr = err / std::max(c, static_cast<double>(exp));
    std::fprintf(stderr, "%d-%g card, expecting %zu, %%%g error\n", isz, double(c), exp, ferr * 100.);
    err = std::abs(ocard2 - exp);
    ferr = err / std::max(ocard2, static_cast<double>(exp));
    std::fprintf(stderr, "%d-after SS conversion, %g card, expecting %zu, %%%g/%g error\n", isz, double(ocard2), exp, ferr * 100., err);
    }}
}
