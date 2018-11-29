/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <random>

#include "../detail/common.hpp"
#include "../next.hpp"

namespace rsm {
namespace detail {

template<>
inline float next(next_value_t<float>, std::mt19937& generator)
{
    return u32_as_float(generator());
}

template<>
inline double next(next_value_t<double>, std::mt19937& generator)
{
    return static_cast<double>(u32_as_float(generator()));
}

template<>
inline float next(next_value_t<float>, std::mt19937_64& generator)
{
    return static_cast<float>(u64_as_double(generator()));
}

template<>
inline double next(next_value_t<double>, std::mt19937_64& generator)
{
    return u64_as_double(generator());
}

} // detail
} // rsm
