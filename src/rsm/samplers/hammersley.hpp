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

template<unsigned int MaxDim, typename FloatType=double>
struct hammersley_sampler
{
    static_assert(MaxDim > 0, "Maximum dimension must be greater than zero");

    explicit hammersley_sampler(size_t max_samples, unsigned int dim=0, uint64_t offset=0)
        : base_dim(dim)
        , offset(offset)
    {
        assert(max_samples > 0);
        inv_max_samples = FloatType(1.0) / max_samples;

        const auto& g_primes = detail::primes_t::get();
        const auto& g_permutations = detail::lds_permutations_t::get();

        assert(dim + MaxDim-1 <= g_primes.N);
        for(unsigned int i=0; i<MaxDim-1; ++i) {
            base[i] = g_primes.p[dim+i];
            permutation[i] = &g_permutations.p[g_primes.sum[dim+i]];
        }
    }

    size_t max_samples() const
    {
        return static_cast<size_t>(FloatType(1.0) / inv_max_samples);
    }

    FloatType inv_max_samples;
    std::array<uint32_t, MaxDim-1> base;
    std::array<const uint16_t*, MaxDim-1> permutation;
    unsigned int base_dim;
    mutable uint64_t offset;
};

namespace detail {

template<typename T, unsigned int MaxDim>
T sample_hammersley(const hammersley_sampler<MaxDim>& sampler, unsigned int dim, uint64_t offset)
{
    assert(dim < MaxDim);
    assert(offset * sampler.inv_max_samples <= 1.0);

    if(dim == 0) {
        return detail::variate<T>(offset * T(sampler.inv_max_samples));
    }
    else {
        unsigned int dim_offset = sampler.base_dim + dim - 1;
        return detail::variate<T>(radical_inverse<T>(dim_offset, sampler.base[dim-1], sampler.permutation[dim-1], offset));
    }
}

} // detail

template<typename T, unsigned int MaxDim>
T sample(const hammersley_sampler<MaxDim>& sampler)
{
    return detail::sample_hammersley<T>(sampler, 0, sampler.offset++);
}

template<unsigned int N, typename T, unsigned int MaxDim>
T sample_vec(const hammersley_sampler<MaxDim>& sampler)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");

    T v;
    using Scalar = typename std::decay<decltype(v[0])>::type;
    for(unsigned int dim=0; dim<N; ++dim) {
        v[dim] = detail::sample_hammersley<Scalar>(sampler, dim, sampler.offset);
    }
    ++sampler.offset;
    return v;
}

template<unsigned int N, typename T, unsigned int MaxDim>
void sample(const hammersley_sampler<MaxDim>& sampler, T* buffer, size_t count=0)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype(buffer[0])>::type;

    size_t requested_samples = (count > 0) ? count : sampler.max_samples();
    assert(sampler.offset + requested_samples <= sampler.max_samples());

    size_t output_index = 0;
    for(size_t i=0; i<requested_samples; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[output_index++] = detail::sample_hammersley<Scalar>(sampler, dim, sampler.offset);
        }
        ++sampler.offset;
    }
}

template<typename T, unsigned int MaxDim>
void sample(const hammersley_sampler<MaxDim>& sampler, T* buffer, size_t count=0)
{
    sample<1>(sampler, buffer, count);
}

template<unsigned int N, typename T, unsigned int MaxDim>
void sample_vec(const hammersley_sampler<MaxDim>& sampler, T* buffer, size_t count=0)
{
    static_assert(N > 0 && N <= MaxDim, "Requested number of dimensions is not in valid range");
    using Scalar = typename std::decay<decltype((*buffer)[0])>::type;

    size_t requested_samples = (count > 0) ? count : sampler.max_samples();
    assert(sampler.offset + requested_samples <= sampler.max_samples());

    size_t output_index = 0;
    for(size_t i=0; i<requested_samples; ++i) {
        for(unsigned int dim=0; dim<N; ++dim) {
            buffer[output_index][dim] = detail::sample_hammersley<Scalar>(sampler, dim, sampler.offset);
        }
        ++sampler.offset;
        ++output_index;
    }
}

} // rsm
