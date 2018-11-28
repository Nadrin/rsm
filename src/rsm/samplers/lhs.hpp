/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <type_traits>

#include "../next.hpp"
#include "../options.hpp"
#include "../utils.hpp"

namespace rsm {

struct lhs_sampler
{
    explicit lhs_sampler(options_t options=opt::jitter)
        : options(options)
    {}
    options_t options;
};

template<unsigned int N, typename T, typename Generator>
void sample(const lhs_sampler& sampler, Generator& generator, T* buffer, size_t count)
{
    static_assert(N > 0, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype(buffer[0])>::type;

    if(count == 0) {
        return;
    }

    size_t output_index = 0;
    Scalar delta = Scalar(1.0) / count;
    if(sampler.options & opt::jitter) {
        for(size_t i=0; i<count; ++i) {
            for(unsigned int dim=0; dim<N; ++dim) {
                Scalar jitter = next<Scalar>(generator);
                buffer[output_index++] = detail::variate<Scalar>((i + jitter) * delta);
            }
        }
    }
    else {
        for(size_t i=0; i<count; ++i) {
            for(unsigned int dim=0; dim<N; ++dim) {
                buffer[output_index++] = detail::variate<Scalar>((i + Scalar(0.5)) * delta);
            }
        }
    }
    shuffle_inner<N>(generator, &buffer[0], &buffer[N * count]);
}

template<typename T, typename Generator>
void sample(const lhs_sampler& sampler, Generator& generator, T* buffer, size_t count)
{
    sample<1>(sampler, generator, buffer, count);
}

template<unsigned int N, typename T, typename Generator>
void sample_vec(const lhs_sampler& sampler, Generator& generator, T* buffer, size_t count)
{
    static_assert(N > 0, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype((*buffer)[0])>::type;

    if(count == 0) {
        return;
    }

    Scalar delta = Scalar(1.0) / count;
    if(sampler.options & opt::jitter) {
        for(size_t i=0; i<count; ++i) {
            for(unsigned int dim=0; dim<N; ++dim) {
                Scalar jitter = next<Scalar>(generator);
                buffer[i][dim] = detail::variate<Scalar>((i + jitter) * delta);
            }
        }
    }
    else {
        for(size_t i=0; i<count; ++i) {
            for(unsigned int dim=0; dim<N; ++dim) {
                buffer[i][dim] = detail::variate<Scalar>((i + Scalar(0.5)) * delta);
            }
        }
    }
    for(unsigned int dim=0; dim<N; ++dim) {
        for(size_t i=0; i<count; ++i) {
            uint32_t r = next(generator, i, count);
            std::swap(buffer[i][dim], buffer[r][dim]);
        }
    }
}

} // rsm
