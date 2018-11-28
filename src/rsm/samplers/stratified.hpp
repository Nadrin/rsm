/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <array>
#include <utility>

#include "../detail/common.hpp"
#include "../next.hpp"
#include "../range.hpp"
#include "../options.hpp"
#include "../utils.hpp"

namespace rsm {

template<unsigned int MaxDim>
struct stratified_sampler
{
    static_assert(MaxDim > 0, "Maximum dimension must be greater than zero");

    explicit stratified_sampler(unsigned int strata_count, options_t options=opt::jitter)
        : options(options)
    {
        strata.fill(strata_count);
    }
    explicit stratified_sampler(const std::array<unsigned int, MaxDim>& strata, options_t options=opt::jitter)
        : strata(strata)
        , options(options)
    {}

    size_t total_strata(unsigned int N=MaxDim) const
    {
        assert(N > 0 && N <= MaxDim);
        unsigned int count = strata[0];
        for(unsigned int i=1; i<N; ++i) {
            count *= strata[i];
        }
        return count;
    }

    std::array<unsigned int, MaxDim> strata;
    options_t options;
};

template<unsigned int N, typename T, unsigned int MaxDim>
range_t<N, T> range(T* buffer, const stratified_sampler<MaxDim>& sampler, uint16_t stride=1)
{
    std::array<size_t, N> range_size;
    for(unsigned int i=0; i<N; ++i) {
        range_size[i] = sampler.strata[i];
    }
    return range_t<N, T>{buffer, range_size, stride};
}

template<unsigned int N, typename T, typename Generator, unsigned int MaxDim>
void sample(const stratified_sampler<MaxDim>& sampler, Generator& generator, T* buffer, size_t count=0)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions out of range");
    using Scalar = typename std::decay<decltype(buffer[0])>::type;

    size_t requested_samples = (count > 0) ? count : sampler.total_strata(N);

    // If # of requested samples matches strata configuration do regular sampling.
    if(requested_samples == sampler.total_strata(N)) {
        Scalar delta[N];
        for(unsigned int dim=0; dim<N; ++dim) {
            delta[dim] = Scalar(1.0) / sampler.strata[dim];
        }
        if(sampler.options & opt::jitter) {
            for(auto it : range<N>(buffer, sampler, N)) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    Scalar jitter = next<Scalar>(generator);
                    it.value[dim] = detail::variate<Scalar>((it.index[dim] + jitter) * delta[dim]);
                }
            }
        }
        else {
            for(auto it : range<N>(buffer, sampler, N)) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    it.value[dim] = detail::variate<Scalar>((it.index[dim] + Scalar(0.5)) * delta[dim]);
                }
            }
        }
        if(sampler.options & opt::shuffle) {
            shuffle<N>(generator, &buffer[0], &buffer[N * requested_samples]);
        }
    }
    // Otherwise perform Latin hypercube sampling (LHS).
    else {
        size_t output_index = 0;
        Scalar delta = Scalar(1.0) / requested_samples;
        if(sampler.options & opt::jitter) {
            for(size_t i=0; i<requested_samples; ++i) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    Scalar jitter = next<Scalar>(generator);
                    buffer[output_index++] = detail::variate<Scalar>((i + jitter) * delta);
                }
            }
        }
        else {
            for(size_t i=0; i<requested_samples; ++i) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    buffer[output_index++] = detail::variate<Scalar>((i + Scalar(0.5)) * delta);
                }
            }
        }
        shuffle_inner<N>(generator, &buffer[0], &buffer[N * requested_samples]);
    }
}

template<typename T, typename Generator, unsigned int MaxDim>
void sample(const stratified_sampler<MaxDim>& sampler, Generator& generator, T* buffer, size_t count=0)
{
    sample<1>(sampler, generator, buffer, count);
}

template<unsigned int N, typename T, typename Generator, unsigned int MaxDim>
void sample_vec(const stratified_sampler<MaxDim>& sampler, Generator& generator, T* buffer, size_t count=0)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions out of range");
    using Scalar = typename std::decay<decltype((*buffer)[0])>::type;

    size_t requested_samples = (count > 0) ? count : sampler.total_strata(N);

    // If # of requested samples matches strata configuration do regular sampling.
    if(requested_samples == sampler.total_strata(N)) {
        Scalar delta[N];
        for(unsigned int dim=0; dim<N; ++dim) {
            delta[dim] = Scalar(1.0) / sampler.strata[dim];
        }
        if(sampler.options & opt::jitter) {
            for(auto it : range<N>(buffer, sampler)) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    Scalar jitter = next<Scalar>(generator);
                    (*it.value)[dim] = detail::variate<Scalar>((it.index[dim] + jitter) * delta[dim]);
                }
            }
        }
        else {
            for(auto it : range<N>(buffer, sampler)) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    (*it.value)[dim] = detail::variate<Scalar>((it.index[dim] + Scalar(0.5)) * delta[dim]);
                }
            }
        }
        if(sampler.options & opt::shuffle) {
            shuffle(generator, &buffer[0], &buffer[requested_samples]);
        }
    }
    // Otherwise perform Latin hypercube sampling (LHS).
    else {
        Scalar delta = Scalar(1.0) / requested_samples;
        if(sampler.options & opt::jitter) {
            for(size_t i=0; i<requested_samples; ++i) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    Scalar jitter = next<Scalar>(generator);
                    buffer[i][dim] = detail::variate<Scalar>((i + jitter) * delta);
                }
            }
        }
        else {
            for(size_t i=0; i<requested_samples; ++i) {
                for(unsigned int dim=0; dim<N; ++dim) {
                    buffer[i][dim] = detail::variate<Scalar>((i + Scalar(0.5)) * delta);
                }
            }
        }
        for(unsigned int dim=0; dim<N; ++dim) {
            for(size_t i=0; i<requested_samples; ++i) {
                uint32_t r = next(generator, i, requested_samples);
                std::swap(buffer[i][dim], buffer[r][dim]);
            }
        }
    }
}

} // rsm
