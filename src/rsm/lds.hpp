/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>

#include "detail/lds.hpp"

namespace rsm {

template<typename T>
T radical_inverse(uint16_t base, uint64_t value)
{
    assert(base >= 2);
    const T inv_base = T(1) / base;
    T inv_base_n = T(1);
    uint64_t inverse = 0;
    for(uint64_t n; value > 0; value = n) {
        n = value / base;
        uint16_t d = static_cast<uint16_t>(value - n * base);
        inverse = inverse * base + d;
        inv_base_n *= inv_base;
    }
    return inverse * inv_base_n;
}

template<typename T>
T radical_inverse_scrambled(uint16_t base, const uint16_t* perm, uint64_t value)
{
    assert(base >= 2);
    const T inv_base = T(1) / base;
    T inv_base_n = T(1);
    uint64_t inverse = 0;
    for(uint64_t n; value > 0; value = n) {
        n = value / base;
        uint16_t d = static_cast<uint16_t>(value - n * base);
        inverse = inverse * base + perm[d];
        inv_base_n *= inv_base;
    }
    return inverse * inv_base_n;
}

} // rsm
