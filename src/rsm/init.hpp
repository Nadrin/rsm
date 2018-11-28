/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>

#include "detail/primes.hpp"
#include "detail/ldsperm.hpp"

#include "generators/pcg32.hpp"

#ifndef RSM_MAX_LDS_DIMENSIONS
#define RSM_MAX_LDS_DIMENSIONS 128
#endif

static_assert(RSM_MAX_LDS_DIMENSIONS >= 2, "RSM_MAX_LDS_DIMENSIONS must be at least 2");

namespace rsm {

inline bool init(const allocator_t& allocator, uint16_t max_lds_dimensions = RSM_MAX_LDS_DIMENSIONS)
{
    assert(max_lds_dimensions >= 2);
    if(!detail::primes_t::get().initialize(max_lds_dimensions, allocator)) {
        return false;
    }

    pcg32 g;
    if(!detail::lds_permutations_t::get().initialize(allocator, g)) {
        return false;
    }
    return true;
}

inline bool init(uint16_t max_lds_dimensions = RSM_MAX_LDS_DIMENSIONS)
{
    return init(detail::default_allocator(), max_lds_dimensions);
}

inline void shutdown(const allocator_t& allocator)
{
    detail::lds_permutations_t::get().free(allocator);
    detail::primes_t::get().free(allocator);
}

inline void shutdown()
{
    shutdown(detail::default_allocator());
}

} // rsm
