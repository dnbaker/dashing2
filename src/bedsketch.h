#pragma once
#ifndef DASHING2_BEDSKETCH_H__
#define DASHING2_BEDSKETCH_H__
#include "d2.h"

namespace dashing2 {

std::vector<RegT> bed2sketch(std::string path, const Dashing2Options &opts);
}

#endif
