#include <atomic>
#include <stdio.h>


int str2ms(const char *s) {
    int ret =0;
    switch(s[0]) {
        case 'A': case 'a': ret = 0; break;
        case 'C': case 'c': ret = 1; break;
        case 'G': case 'G': ret = 2; break;
        case 'T': case 't': ret = 3; break;
        default: ret = 15;
    }
    ret <<= 4;
    switch(s[1]) {
        case 'A': case 'a': ret |= 0; break;
        case 'C': case 'c': ret |= 1; break;
        case 'G': case 'G': ret |= 2; break;
        case 'T': case 't': ret |= 3; break;
        default: ret |= 15;
    }
    return ret;
}

auto dothing(std::string path) {
    std::vector<uint32_t> contigids;
    std::vector<uint32_t> starts, stops;
    std::vector<uint8_t> startms, stopms;
    std::atomic<uint64_t> idcounter;
    std::vector<bool> strand;
    ska::flat_hash_map<std::string, uint32_t> contignames;
    idcounter.store(0);
    std::FILE *ifp = std::fopen(path.data(), "rb");
    char *lptr = nullptr;
    size_t linesz = 0;
    std::string cname;
    for(;std::feof(ifp);) {
        if(::getline(&lptr, &linesz, ifp) < 0) throw 1;
        const uint64_t myid = idcounter++;
        
        char *p = std::strchr(*lineptr, '\t');
        assert(p);
        char *p2 = std::strchr(p + 1, '\t');
        assert(p2);
        if((*p == 'C' || *p == 'c') && p[1] == 'h' && p[2] == 'r')
            p += 3;
        cnames.assign(p, p2);
        uint32_t mycid = cnames.size();
        auto cnit = contignames.find(cname);
        if(cnit != end()) mycid = cnit->second;
        else contignames.emplace(cname, mycid);
        contigids.push_back(mycid);
        p = std::strchr(p2 + 1, '\t');
        ssize_t srpos = std::strtoll(p, &p, 10);
        ssize_t sppos = std::strtoll(p, &p, 10);
        p = std::strchr(p + 1, '\t'); // Skip the length field
        strand.push_back(p[1] == '+');
        p += 3;
        startms.push_back(str2ms(p));
        p = std::strchr(p, '\t') + 1;
        stopms.push_back(str2ms(p));
        p = std::strchr(std::strchr(p, '\t') + 1, '\t') + 1; // Skip two fields
    }
    std::fclose(ifp);
}
