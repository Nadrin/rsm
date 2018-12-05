/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include "../detail/common.hpp"

namespace rsm {

template<unsigned int N, typename U>
void ncube(const U* min, const U* max, const U* u, U* p)
{
    static_assert(N > 0, "Number of dimensions must be greater than zero");

    for(unsigned int i=0; i<N; ++i) {
        p[i] = min[i] + u[i] * (max[i] - min[i]);
    }
}

template<unsigned int N, typename U>
U ncube_pdf(const U* min, const U* max)
{
    static_assert(N > 0, "Number of dimensions must be greater than zero");

    U volume = U(1.0);
    for(unsigned int i=0; i<N; ++i) {
        volume *= (max[i] - min[i]);
    }
    return U(1.0) / volume;
}

template<unsigned int N, typename T>
T ncube(const T& min, const T& max, const T& u)
{
    static_assert(N > 0, "Number of dimensions must be greater than zero");

    T p;
    for(unsigned int i=0; i<N; ++i) {
        p[i] = min[i] + u[i] * (max[i] - min[i]);
    }
    return p;
}

template<unsigned int N, typename T, typename U=typename detail::component<T>::type>
U ncube_pdf(const T& min, const T& max)
{
    static_assert(N > 0, "Number of dimensions must be greater than zero");

    U volume = U(1.0);
    for(unsigned int i=0; i<N; ++i) {
        volume *= (max[i] - min[i]);
    }
    return U(1.0) / volume;
}

template<typename U>
void rectangle(const U* min, const U* max, const U* u, U* p)
{
    ncube<2>(min, max, u, p);
}

template<typename U>
U rectangle_pdf(const U* min, const U* max)
{
    return U(1.0) / ((max[0] - min[0]) * (max[1] - min[1]));
}

template<typename T>
T rectangle(const T& min, const T& max, const T& u)
{
    return ncube<2>(min, max, u);
}

template<typename T, typename U=typename detail::component<T>::type>
U rectangle_pdf(const T& min, const T& max)
{
    return U(1.0) / ((max[0] - min[0]) * (max[1] - min[1]));
}

template<typename T, typename U=typename detail::component<T>::type>
T rectangle(const T& min, const T& max, U u1, U u2)
{
    T p;
    p[0] = min[0] + u1 * (max[0] - min[0]);
    p[1] = min[1] + u2 * (max[1] - min[1]);
    return p;
}

} // rsm
