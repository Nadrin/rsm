/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <limits>
#include <algorithm>

namespace rsm {
namespace detail {

template<typename T> constexpr T one_minus_eps() = delete;
template<> constexpr float  one_minus_eps() { return 0.99999994f; }
template<> constexpr double one_minus_eps() { return 0.99999999999999989; }

template<typename T, typename U> constexpr T inv_max()
{
    return 1.0 / std::numeric_limits<U>::max();
}

template<typename T>
T variate(T value)
{
    return std::min<T>(value, one_minus_eps<T>());
}

} // detail
} // rsm
