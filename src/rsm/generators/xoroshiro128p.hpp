/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

/*
 * This file contains an implementation of xoroshiro128+ 1.0 random number generator adapted for RSM.
 * Original C implementation written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * and released into public domain.
 * Re-licensed under the MIT license for use within RSM library.
 */

/*
 * This is xoroshiro128+ 1.0, our best and fastest small-state generator
 * for floating-point numbers. We suggest to use its upper bits for
 * floating-point generation, as it is slightly faster than
 * xoroshiro128**. It passes all tests we are aware of except for the four
 * lower bits, which might fail linearity tests (and just those), so if
 * low linear complexity is not considered an issue (as it is usually the
 * case) it can be used to generate 64-bit outputs, too; moreover, this
 * generator has a very mild Hamming-weight dependency making our test
 * (http://prng.di.unimi.it/hwd.php) fail after 8 TB of output; we believe
 * this slight bias cannot affect any application. If you are concerned,
 * use xoroshiro128** or xoshiro256+.
 * NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
 * better results in our test than the 2016 version (a=55, b=14, c=36).
 */

#pragma once

#include <cstdint>
#include <limits>

#include "../detail/common.hpp"
#include "../next.hpp"

#include "splitmix64.hpp"

namespace rsm {

class xoroshiro128p
{
public:
    using result_type = uint64_t;

    xoroshiro128p()
        : m_state{0x6ca496188ffdcf87ull, 0x492492d4cd7e9082ull}
    {}

    xoroshiro128p(uint64_t s)
    {
        seed(s);
    }

    void seed(uint64_t s)
    {
        splitmix64 g(s);
        m_state[0] = g();
        m_state[1] = g();
    }

    result_type operator()()
    {
        uint64_t s0 = m_state[0];
        uint64_t s1 = m_state[1];
        uint64_t result = s0 + s1;

        s1 ^= s0;
        m_state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        m_state[1] = rotl(s1, 37); // c

        return result;
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
    static inline uint64_t rotl(uint64_t x, int k)
    {
        return (x << k) | (x >> (64 - k));
    }
    uint64_t m_state[2];
};

namespace detail {

template<>
inline float next(next_value_t<float>, xoroshiro128p& generator)
{
    return u32_as_float(static_cast<uint32_t>(generator()));
}

template<>
inline double next(next_value_t<double>, xoroshiro128p& generator)
{
    return u64_as_double(generator());
}

} // detail
} // rsm
