/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

/*
 * This file contains an implementation of xoroshiro64* 1.0 random number generator adapted for RSM.
 * Original C implementation written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * and released into public domain.
 * Re-licensed under the MIT license for use within RSM library.
 */

/*
 * This is xoroshiro64* 1.0, our best and fastest 32-bit small-state generator
 * for 32-bit floating-point numbers. We suggest to use its upper bits for
 * floating-point generation, as it is slightly faster than
 * xoroshiro64**. It passes all tests we are aware of except for
 * linearity tests, as the lowest six bits have low linear complexity, so
 * if low linear complexity is not considered an issue (as it is usually
 * the case) it can be used to generate 32-bit outputs, too.
 */

#pragma once

#include <cstdint>
#include <limits>

#include "../detail/common.hpp"
#include "../next.hpp"

#include "splitmix64.hpp"

namespace rsm {

class xoroshiro64s
{
public:
    using result_type = uint32_t;

    xoroshiro64s()
        : m_state{0xcb308ebeul, 0xd97012dcul}
    {}

    xoroshiro64s(uint64_t s)
    {
        seed(s);
    }

    void seed(uint64_t s)
    {
        uint64_t state64 = splitmix64(s)();
        m_state[0] = static_cast<uint32_t>(state64);
        m_state[1] = static_cast<uint32_t>(state64 >> 32);
    }

    result_type operator()()
    {
        uint32_t s0 = m_state[0];
        uint32_t s1 = m_state[1];
        uint32_t result_star = s0 * 0x9e3779bbul;

        s1 ^= s0;
        m_state[0] = rotl(s0, 26) ^ s1 ^ (s1 << 9); // a, b
        m_state[1] = rotl(s1, 13); // c

        return result_star;
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
    static inline uint32_t rotl(uint32_t x, int k)
    {
        return (x << k) | (x >> (32 - k));
    }
    uint32_t m_state[2];
};

namespace detail {

template<>
inline float next(next_value_t<float>, xoroshiro64s& generator)
{
    return u32_as_float(generator());
}

template<>
inline double next(next_value_t<double>, xoroshiro64s& generator)
{
    return static_cast<double>(u32_as_float(generator()));
}

} // detail
} // rsm
