/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

/*
 * This file contains an implementation of SplitMix64 random number generator adapted for RSM.
 * Original C implementation written in 2015 by Sebastiano Vigna (vigna@acm.org)
 * and released into public domain.
 * Re-licensed under the MIT license for use within RSM library.
 */

/*
 * This is a fixed-increment version of Java 8's SplittableRandom generator
 * See http://dx.doi.org/10.1145/2714064.2660195 and
 * http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 *
 * It is a very fast generator passing BigCrush, and it can be useful if
 * for some reason you absolutely want 64 bits of state; otherwise, we
 * rather suggest to use a xoroshiro128+ (for moderately parallel
 * computations) or xorshift1024* (for massively parallel computations)
 * generator.
 */

#pragma once

#include <cstdint>
#include <limits>

#include "../detail/common.hpp"
#include "../next.hpp"

namespace rsm {

class splitmix64
{
public:
    using result_type = uint64_t;

    splitmix64()
        : m_state(0x27c6003152ca78dull)
    {}

    splitmix64(uint64_t state)
        : m_state(state)
    {}

    void seed(uint64_t state)
    {
        m_state = state;
    }

    result_type operator()()
    {
        uint64_t z = (m_state += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }

    static constexpr result_type min()
    {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type max()
    {
        return std::numeric_limits<result_type>::max();
    }

private:
    uint64_t m_state;
};

namespace detail {

template<>
inline float next(next_value_t<float>, splitmix64& generator)
{
    return u32_as_float(static_cast<uint32_t>(generator()));
}

template<>
inline double next(next_value_t<double>, splitmix64& generator)
{
    return u64_as_double(generator());
}

} // detail
} // rsm
