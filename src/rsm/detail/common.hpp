/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>
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

// Make sure random variate is kept in [0..1) range.
template<typename T>
T variate(T value)
{
    return std::min<T>(value, one_minus_eps<T>());
}

// Fast cast of 32-bit unsigned int to float in [0..1) interval.
// See: http://xoshiro.di.unimi.it/#remarks, "Generating uniform doubles in the unit interval"
inline float u32_as_float(uint32_t value)
{
    union {
        uint32_t u;
        float f;
    } r;
    r.u = 0x3f800000ul | (value >> 9);
    return r.f - 1.0f;
}

// Fast cast of 64-bit unsigned int to double in [0..1) interval.
// See: http://xoshiro.di.unimi.it/#remarks, "Generating uniform doubles in the unit interval"
inline double u64_as_double(uint64_t value)
{
    union {
        uint64_t u;
        double d;
    } r;
    r.u = 0x3ff0000000000000ull | (value >> 12);
    return r.d - 1.0;
}

} // detail
} // rsm
