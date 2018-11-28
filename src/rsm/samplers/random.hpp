/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <type_traits>

#include "../next.hpp"

namespace rsm {

template<typename T, typename Generator>
T sample(Generator& generator)
{
    return next<T>(generator);
}

template<unsigned int N, typename T, typename Generator>
T sample_vec(Generator& generator)
{
    static_assert(N > 0, "Requested number of dimensions is not in valid range");

    T v;
    using Scalar = typename std::decay<decltype(v[0])>::type;
    for(unsigned int dim=0; dim<N; ++dim) {
        v[dim] = next<Scalar>(generator);
    }
    return v;
}

template<unsigned int N, typename T, typename Generator>
void sample(Generator& generator, T* buffer, size_t count)
{
    static_assert(N > 0, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype(buffer[0])>::type;

    size_t output_index = 0;
    for(size_t i=0; i<count; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[output_index++] = next<Scalar>(generator);
        }
    }
}

template<typename T, typename Generator>
void sample(Generator& generator, T* buffer, size_t count)
{
    sample<1>(generator, buffer, count);
}

template<unsigned int N, typename T, typename Generator>
void sample_vec(Generator& generator, T* buffer, size_t count)
{
    static_assert(N > 0, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype((*buffer)[0])>::type;

    for(size_t i=0; i<count; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[i][dim] = next<Scalar>(generator);
        }
    }
}

} // rsm
