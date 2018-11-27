/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>

#include "detail/common.hpp"

namespace rsm {

namespace detail {

template<typename T> struct next_value_t{};

template<typename Generator>
inline uint32_t next(next_value_t<uint32_t>, Generator& generator)
{
    return static_cast<uint32_t>(generator());
}

template<typename Generator>
inline uint32_t next(next_value_t<uint32_t>, Generator& generator, uint32_t range)
{
    // Based on http://www.pcg-random.org/posts/bounded-rands.html
    // Lemire's method with t-opt.
    uint32_t x = static_cast<uint32_t>(generator() - Generator::min());
    uint64_t m = uint64_t(x) * uint64_t(range);
    uint32_t l = uint32_t(m);
    if(l < range) {
        uint32_t t = (-range) % range;
        while(l < t) {
            x = static_cast<uint32_t>(generator() - Generator::min());
            m = uint64_t(x) * uint64_t(range);
            l = uint32_t(m);
        }
    }
    return m >> 32;
}

template<typename Generator>
inline uint64_t next(next_value_t<uint64_t>, Generator& generator)
{
    return static_cast<uint64_t>(generator());
}

template<typename Generator>
inline uint64_t next(next_value_t<uint64_t>, Generator& generator, uint64_t range)
{
    // Based on http://www.pcg-random.org/posts/bounded-rands.html
    // Debiased modulo (once) method.
    // Lemire's method is significantly faster but in 64-bit output case it uses 128-bit integeres.
    // We can't rely on those being available (especially when compiling for CUDA).
    uint64_t x, r;
    do {
        x = static_cast<uint64_t>(generator() - Generator::min());
        r = x % range;
    } while(x - r > uint64_t(-range));
    return r;
}

template<typename Generator>
inline float next(next_value_t<float>, Generator& generator)
{
    // This method introduces slight bias but we need it to be as fast as possible.
    // See: http://mumble.net/~campbell/tmp/random_real.c
    constexpr float inv_range = 1.0 / (Generator::max() - Generator::min());
    return detail::variate<float>((generator() - Generator::min()) * inv_range);
}

template<typename Generator>
inline double next(next_value_t<double>, Generator& generator)
{
    // This method introduces slight bias but we need it to be as fast as possible.
    // See: http://mumble.net/~campbell/tmp/random_real.c
    constexpr double inv_range = 1.0 / (Generator::max() - Generator::min());
    return detail::variate<double>((generator() - Generator::min()) * inv_range);
}

} // detail

template<typename T, typename Generator>
inline T next(Generator& generator)
{
    return detail::next<Generator>(detail::next_value_t<T>{}, generator);
}

template<typename T, typename Generator>
inline T next(Generator& generator, T min, T max)
{
    assert(min < max);
    return min + detail::next<Generator>(detail::next_value_t<T>{}, generator, max - min);
}

} // rsm
