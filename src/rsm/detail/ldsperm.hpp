/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>

#include "common.hpp"
#include "primes.hpp"
#include "../utils.hpp"

namespace rsm {
namespace detail {

struct lds_permutations_t
{
    uint32_t N = 0;
    uint16_t* p = nullptr;

    template<typename Generator>
    bool initialize(const allocator_t& allocator, Generator& g)
    {
        if(N > 0) {
            assert(p);
            return true;
        }

        const auto& primes = primes_t::get();
        N = primes.sum[primes.N];
        p = detail::alloc<uint16_t>(allocator, N);
        if(!p) {
            return false;
        }

        uint16_t* p_offset = p;
        for(uint32_t i=0; i<primes.N; ++i) {
            const uint32_t pn = primes.p[i];
            for(uint32_t j=0; j<pn; ++j) {
                p_offset[j] = static_cast<uint16_t>(j);
            }
            // Don't shuffle first two dimensions as the distribution is good even without scrambling.
            if(i >= 2) {
                // Shuffle to generate a scrambling permutation (but always map 0 to 0).
                shuffle(g, &p_offset[1], &p_offset[pn]);
            }
            p_offset += pn;
        }
        return true;
    }

    void free(const allocator_t& allocator)
    {
        N = 0;
        detail::free(allocator, p);
    }

    static lds_permutations_t& get()
    {
        static lds_permutations_t self;
        return self;
    }
};

} // detail
} // rsm
