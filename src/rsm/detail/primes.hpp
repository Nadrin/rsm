/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <cmath>

#include "memory.hpp"

namespace rsm {
namespace detail {

struct primes_t
{
    uint16_t N = 0;
    uint32_t* p = nullptr;
    uint32_t* sum = nullptr;

    bool initialize(uint16_t num, const allocator_t& allocator)
    {
        if(N > 0) {
            assert(p && sum);
            return true;
        }

        N = num;
        p = detail::alloc<uint32_t>(allocator, N);
        if(!p) {
            return false;
        }
        sum = detail::alloc<uint32_t>(allocator, N+1);
        if(!sum) {
            detail::free(allocator, p);
            return false;
        }

        p[0] = 2;
        p[1] = 3;
        uint16_t i=2;
        for(uint32_t np = p[i-1] + 2; i<N; np += 2) {
            bool is_prime = true;
            for(unsigned int divisor = 2; divisor <= static_cast<unsigned int>(std::sqrt(np)); ++divisor) {
                if((np % divisor) == 0) {
                    is_prime = false;
                    break;
                }
            }
            if(is_prime) {
                p[i++] = np;
            }
        }

        sum[0] = 0;
        for(unsigned int i=1; i<=N; ++i) {
            sum[i] = sum[i-1] + p[i-1];
        }
        return true;
    }

    void free(const allocator_t& allocator)
    {
        N = 0;
        detail::free(allocator, p);
        detail::free(allocator, sum);
    }

    static primes_t& get()
    {
        static primes_t self;
        return self;
    }
};

} // detail
} // rsm
