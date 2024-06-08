#pragma once
#ifndef FASTXMERGE_H__
#define FASTXMERGE_H__

#include "fastxsketch.h"

namespace dashing2 {

class SketchingResult;

SketchingResult merge(SketchingResult *start, size_t n, const std::vector<std::string> &names = std::vector<std::string>());

std::string makedest(Dashing2Options &opts, const std::string &path, bool iskmer=false);

}

#endif // FASTXMERGE_H__
