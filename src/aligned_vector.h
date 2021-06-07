#pragma once

namespace aligned {

template<typename T>
struct vector: public std::vector<T, sse::AlignedAllocator<T, sse::Alignment::AVX512>> {};

}
