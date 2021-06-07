#pragma once
#ifndef DASHING2_LFSKETCH_H__
#define DASHING2_LFSKETCH_H__
#include "d2.h"

namespace dashing2 {

struct LFResult: public std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<RegT>, std::vector<std::string>, std::vector<size_t>>
{
    auto &splice_sites() {return std::get<0>(*this);}
    const auto &splice_sites() const {return std::get<0>(*this);}
    auto &sample_names() {return std::get<1>(*this);}
    const auto &sample_names() const {return std::get<1>(*this);}
    auto &registers() {return std::get<2>(*this);}
    const auto &registers() const {return std::get<2>(*this);}
    auto &filenames() {return std::get<3>(*this);}
    const auto &filenames() const {return std::get<3>(*this);}
    auto &nsamples_per_file() {return std::get<4>(*this);}
    const auto &nsamples_per_file() const {return std::get<4>(*this);}
    static LFResult merge_results(const LFResult *start, size_t n);

};

LFResult lf2sketch(std::string path, const Dashing2Options &opts);
LFResult lf2sketch(std::vector<std::string> paths, const Dashing2Options &opts);
}

#endif
