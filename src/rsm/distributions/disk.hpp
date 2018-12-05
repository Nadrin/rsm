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
void disk(U radius, U u1, U u2, U& px, U& py)
{
    U r = radius * std::sqrt(u1);
    U theta = U(2.0) * detail::pi<U>() * u2;
    px = r * std::cos(theta);
    py = r * std::sin(theta);
}

template<typename U>
constexpr U disk_pdf(U radius)
{
    U inv_radius_sqr = U(1.0) / (radius * radius);
    return detail::inv_pi<U>() * inv_radius_sqr;
}

template<typename T, typename U=typename detail::component<T>::type>
T disk(U radius, U u1, U u2)
{
    T r;
    disk(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename U=typename detail::component<T>::type>
T disk(U radius, const U& u)
{
    T r;
    disk(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void disk(U radius, const U* u, U* p)
{
    disk(radius, u[0], u[1], p[0], p[1], p[2]);
}

template<typename U>
void disk_concentric(U radius, U u1, U u2, U& px, U& py)
{
    U remap_u1 = U(2.0) * u1 - U(1.0);
    U remap_u2 = U(2.0) * u2 - U(1.0);

    if(remap_u1 == U(0.0) && remap_u2 == U(0.0)) {
        px = U(0.0);
        py = U(0.0);
    }
    else {
        U r;
        U theta;
        if(std::abs(remap_u1) > std::abs(remap_u2)) {
            r = radius * remap_u1;
            theta = U(0.25) * detail::pi<U>() * (remap_u2 / remap_u1);
        }
        else {
            r = radius * remap_u2;
            theta = U(0.5) * detail::pi<U>() - U(0.25) * detail::pi<U>() * (remap_u1 / remap_u2);
        }
        px = r * std::cos(theta);
        py = r * std::sin(theta);
    }
}

template<typename U>
constexpr U disk_concentric_pdf(U radius)
{
    return disk_pdf(radius);
}

template<typename T, typename U=typename detail::component<T>::type>
T disk_concentric(U radius, U u1, U u2)
{
    T r;
    disk_concentric(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename U=typename detail::component<T>::type>
T disk_concentric(U radius, const U& u)
{
    T r;
    disk_concentric(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void disk_concentric(U radius, const U* u, U* p)
{
    disk_concentric(radius, u[0], u[1], p[0], p[1], p[2]);
}

} // rsm
