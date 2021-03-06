/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

/*
 * This file contains an implementation of PCG32 random number generator adapted for RSM.
 * Original C implementation Copyright (c) 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

#pragma once

#include <cstdint>
#include <limits>

#include "../detail/common.hpp"
#include "../next.hpp"

#include "splitmix64.hpp"

namespace rsm {

class pcg32
{
public:
    using result_type = uint32_t;

    pcg32()
        : m_state(0x853c49e6748fea9bull)
        , m_inc(0xda3e39cb94b95bdbull)
    {}

    pcg32(uint64_t s)
    {
        seed(s);
    }

    pcg32(uint64_t initstate, uint64_t initstream)
    {
        seed(initstate, initstream);
    }

    void seed(uint64_t s)
    {
        splitmix64 g(s);
        uint64_t initstate = g();
        uint64_t initstream = g();
        seed(initstate, initstream);
    }

    void seed(uint64_t initstate, uint64_t initstream)
    {
        m_state = 0;
        m_inc = (initstream << 1u) | 1u;
        (*this)();
        m_state += initstate;
        (*this)();
    }

    result_type operator()()
    {
        uint64_t oldstate = m_state;
        m_state = oldstate * 6364136223846793005ull + m_inc;
        uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
        uint32_t rot = oldstate >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
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
    uint64_t m_inc;
};

namespace detail {

template<>
inline float next(next_value_t<float>, pcg32& generator)
{
    return u32_as_float(generator());
}

template<>
inline double next(next_value_t<double>, pcg32& generator)
{
    return static_cast<double>(u32_as_float(generator()));
}

} // detail
} // rsm
