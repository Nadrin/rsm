/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cmath>

#include "../detail/common.hpp"

namespace rsm {

template<typename U>
void sphere(U radius, U u1, U u2, U& px, U& py, U& pz)
{
    U z   = U(1.0) - U(2.0) * u1;
    U phi = U(2.0) * detail::pi<U>() * u2;
    U sin_theta_r = radius * std::sqrt(U(1.0) - z * z);
    px = std::cos(phi) * sin_theta_r;
    py = std::sin(phi) * sin_theta_r;
    pz = radius * z;
}

template<typename U>
constexpr U sphere_pdf(U radius)
{
    U inv_radius_sqr = U(1.0) / (radius * radius);
    return U(0.25) * detail::inv_pi<U>() * inv_radius_sqr;
}

template<typename T, typename U=typename detail::component<T>::type>
T sphere(U radius, U u1, U u2)
{
    T r;
    sphere(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename U=typename detail::component<T>::type>
T sphere(U radius, const T& u)
{
    T r;
    sphere(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void sphere(U radius, const U* u, U* p)
{
    sphere(radius, u[0], u[1], p[0], p[1], p[2]);
}

} // rsm
