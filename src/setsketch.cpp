#include <setsketch.h>

namespace sketch {
namespace setsketch {
namespace detail {

std::pair<long double, long double> optimal_parameters(const long double maxreg, const long double minreg, const long double q) noexcept {
    const long double b = std::exp(std::log(maxreg / minreg) / q);
    return {b, maxreg / b};
}
}}}
