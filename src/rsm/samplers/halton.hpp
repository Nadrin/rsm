/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <array>
#include <type_traits>

#include "../detail/common.hpp"
#include "../lds.hpp"

namespace rsm {

template<unsigned int MaxDim>
struct halton_sampler
{
    static_assert(MaxDim > 0, "Maximum dimension must be greater than zero");

    explicit halton_sampler(unsigned int dim=0, uint64_t offset=0)
        : base_dim(dim)
        , offset(offset)
    {
        const auto& g_primes = detail::primes_t::get();
        const auto& g_permutations = detail::lds_permutations_t::get();

        assert(dim + MaxDim <= g_primes.N);
        for(unsigned int i=0; i<MaxDim; ++i) {
            base[i] = g_primes.p[dim+i];
            permutation[i] = &g_permutations.p[g_primes.sum[dim+i]];
        }
    }

    std::array<uint32_t, MaxDim> base;
    std::array<const uint16_t*, MaxDim> permutation;
    unsigned int base_dim;
    mutable uint64_t offset;
};

namespace detail {

template<typename T, unsigned int MaxDim>
T sample_halton(const halton_sampler<MaxDim>& sampler, unsigned int dim, uint64_t offset)
{
    assert(dim < MaxDim);
    unsigned int dim_offset = sampler.base_dim + dim;
    return detail::variate<T>(radical_inverse<T>(dim_offset, sampler.base[dim], sampler.permutation[dim], offset));
}

} // detail

template<typename T, unsigned int MaxDim>
T sample(const halton_sampler<MaxDim>& sampler)
{
    return detail::sample_halton<T>(sampler, 0, sampler.offset++);
}

template<unsigned int N, typename T, unsigned int MaxDim>
T sample_vec(const halton_sampler<MaxDim>& sampler)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");

    T v;
    using Scalar = typename std::decay<decltype(v[0])>::type;
    for(unsigned int dim=0; dim<N; ++dim) {
        v[dim] = detail::sample_halton<Scalar>(sampler, dim, sampler.offset);
    }
    ++sampler.offset;
    return v;
}

template<unsigned int N, typename T, unsigned int MaxDim>
void sample(const halton_sampler<MaxDim>& sampler, T* buffer, size_t count)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype(buffer[0])>::type;

    size_t output_index = 0;
    for(size_t i=0; i<count; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[output_index++] = detail::sample_halton<Scalar>(sampler, dim, sampler.offset);
        }
        ++sampler.offset;
    }
}

template<typename T, unsigned int MaxDim>
void sample(const halton_sampler<MaxDim>& sampler, T* buffer, size_t count)
{
    sample<1>(sampler, buffer, count);
}

template<unsigned int N, typename T, unsigned int MaxDim>
void sample_vec(const halton_sampler<MaxDim>& sampler, T* buffer, size_t count)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype((*buffer)[0])>::type;

    size_t output_index = 0;
    for(size_t i=0; i<count; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[output_index][dim] = detail::sample_halton<Scalar>(sampler, dim, sampler.offset);
        }
        ++sampler.offset;
        ++output_index;
    }
}

} // rsm
