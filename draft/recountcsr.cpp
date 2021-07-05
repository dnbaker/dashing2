#include <atomic>
#include <vector>
#include <cassert>
#include <cstring>
#include <tuple>
#include <string>
#include "flat_hash_map/flat_hash_map.hpp"
#include <stdio.h>


int str2ms(const char *s) {
    int ret =0;
    switch(s[0]) {
        case 'A': case 'a': ret = 0; break;
        case 'C': case 'c': ret = 1; break;
        case 'G': case 'g': ret = 2; break;
        case 'T': case 't': ret = 3; break;
        default: ret = 15;
    }
    ret <<= 4;
    switch(s[1]) {
        case 'A': case 'a': ret |= 0; break;
        case 'C': case 'c': ret |= 1; break;
        case 'G': case 'g': ret |= 2; break;
        case 'T': case 't': ret |= 3; break;
        default: ret |= 15;
    }
    return ret;
}

auto parse_file(std::FILE *ifp) {
    std::fprintf(stderr, "Starting\n");
    std::vector<uint32_t> contigids;
    std::vector<uint64_t> ids;
    std::vector<uint32_t> starts, stops;
    std::vector<uint8_t> startms, stopms;
    std::atomic<uint64_t> idcounter;
    std::vector<uint16_t> counts;
    std::vector<uint64_t> indptr{0};
    ska::flat_hash_map<std::string, uint32_t> contignames;
    idcounter.store(0);
    char *lptr = nullptr;
    size_t linesz = 0;
    std::string cname;
    size_t ln = 0;
    ssize_t rc;
	for(ssize_t rc; (rc = ::getline(&lptr, &linesz, ifp)) >= 0;++ln) {
        if(ln % 1 == 65536) std::fprintf(stderr, "Processed %zu lines, last rc is %zd\n", ln, rc);
        const uint64_t myid = idcounter++;
        
        char *p = std::strchr(lptr, '\t');
        assert(p);
        char *p2 = std::strchr(p + 1, '\t');
        assert(p2);
        if((*p == 'C' || *p == 'c') && p[1] == 'h' && p[2] == 'r')
            p += 3;
        cname.assign(p, p2);
        uint32_t mycid = cname.size();
        auto cnit = contignames.find(cname);
        if(cnit != contignames.end()) mycid = cnit->second;
        else contignames.emplace(cname, mycid);
        contigids.push_back(mycid);
        p = p2 + 1;
        ssize_t srpos = std::strtoll(p, &p, 10);
        ssize_t sppos = std::strtoll(p, &p, 10);
        p = std::strchr(p + 1, '\t'); // Skip the length field
        assert(p[1] == '+' || p[1] == '-' || p[1] == '?' || !std::fprintf(stderr, "p: %s, %d\n", p, p[1]));
        //p += 3;
        p = std::strchr(p + 3, '\t') + 1;
        startms.push_back(str2ms(p));
        p = std::strchr(p, '\t') + 1;
        stopms.push_back(str2ms(p));
        p = std::strchr(std::strchr(std::strchr(p, '\t') + 1, '\t') + 1, '\t') + 1; // Skip two fields
        uint64_t id;
        int32_t ct;
        size_t nids = 0;
        for(char *p3 = p;*p3 && *p3 != '\t';) {
            id = std::strtoull(p3 + 1, &p3, 10);
            ct = std::strtoll(p3 + 1, &p3, 10);
            counts.push_back(ct);
            ids.push_back(id);
            ++nids;
        }
        indptr.push_back(nids);
    }
    return std::make_tuple(contigids, contignames, counts, ids, indptr);
}

int main(int argc, char **argv) {
    std::FILE *fp = std::fopen(argc == 1 ? "/dev/stdin": argv[1], "r");
    auto [cids, cnames, counts, ids, indptr] = parse_file(fp);
    std::fclose(fp);

    ska::flat_hash_map<uint64_t, uint16_t> mapper;
    for(const auto id: ids) {
        if(mapper.find(id) == mapper.end()) {
            auto oldsz = mapper.size();
            mapper.emplace(id, oldsz);
        }
    }
    assert(mapper.size() < 65536);
    std::transform(ids.begin(), ids.end(), ids.begin(), [&mapper](auto x) {return mapper[x];});
    fp = std::fopen("parsed.cts.u16", "wb");
    std::fwrite(counts.data(), 2, counts.size(), fp);
    std::fclose(fp);
    fp = std::fopen("parsed.ids.u16", "wb");
    for(const auto id: ids) {
        uint16_t sid(id);
        std::fwrite(&sid, 2, 1, fp);
    }
    std::fclose(fp);
    fp = std::fopen("parsed.indptr.u16", "w");
    std::fwrite(indptr.data(), sizeof(indptr[0]), indptr.size(), fp);
    std::fclose(fp);

    for(const auto &pair: mapper)
        std::fprintf(fp, "%zu:%zu\n", size_t(pair.first), size_t(pair.second));
    std::fclose(fp);

    fp = std::fopen("parsed.remap", "w");
    for(const auto &pair: mapper)
        std::fprintf(fp, "%zu:%zu\n", size_t(pair.first), size_t(pair.second));
    std::fclose(fp);

    return 0;
}
