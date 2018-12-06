/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cmath>

#include "../detail/common.hpp"
#include "disk.hpp"

namespace rsm {

template<typename U>
void hemisphere(U radius, U u1, U u2, U& px, U& py, U& pz)
{
    U phi = U(2.0) * detail::pi<U>() * u2;
    U sin_theta_r = radius * std::sqrt(U(1.0) - u1 * u1);
    px = std::cos(phi) * sin_theta_r;
    py = std::sin(phi) * sin_theta_r;
    pz = radius * u1;
}

template<typename U>
constexpr U hemisphere_pdf(U radius)
{
    U inv_radius_sqr = U(1.0) / (radius * radius);
    return U(0.5) * detail::inv_pi<U>() * inv_radius_sqr;
}

template<typename T, typename U=typename detail::component<T>::type>
T hemisphere(U radius, U u1, U u2)
{
    T r;
    hemisphere(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename V, typename U=typename detail::component<T>::type>
T hemisphere(U radius, const V& u)
{
    T r;
    hemisphere(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void hemisphere(U radius, const U* u, U* p)
{
    hemisphere(radius, u[0], u[1], p[0], p[1], p[2]);
}

template<typename U>
void hemisphere_cosine(U radius, U u1, U u2, U& px, U& py, U& pz)
{
    U disk_px, disk_py;
    disk(U(1.0), u1, u2, disk_px, disk_py);
    px = radius * disk_px;
    py = radius * disk_py;
    pz = radius * std::sqrt(U(1.0) - disk_px*disk_px - disk_py*disk_py);
}

template<typename U>
constexpr U hemisphere_cosine_pdf(U radius, U cos_theta)
{
    U inv_radius_sqr = U(1.0) / (radius * radius);
    return cos_theta * detail::inv_pi<U>() * inv_radius_sqr;
}

template<typename U>
constexpr U hemisphere_cosine_pdf_s(U radius, U cos_theta)
{
    U inv_radius_sqr = U(1.0) / (radius * radius);
    U sin_theta = std::sqrt(U(1.0) - cos_theta);
    return sin_theta * cos_theta * detail::inv_pi<U>() * inv_radius_sqr;
}

template<typename T, typename U=typename detail::component<T>::type>
T hemisphere_cosine(U radius, U u1, U u2)
{
    T r;
    hemisphere_cosine(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename V, typename U=typename detail::component<T>::type>
T hemisphere_cosine(U radius, const V& u)
{
    T r;
    hemisphere_cosine(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void hemisphere_cosine(U radius, const U* u, U* p)
{
    hemisphere_cosine(radius, u[0], u[1], p[0], p[1], p[2]);
}

template<typename U>
void hemisphere_cosine_concentric(U radius, U u1, U u2, U& px, U& py, U& pz)
{
    U disk_px, disk_py;
    disk_concentric(U(1.0), u1, u2, disk_px, disk_py);
    px = radius * disk_px;
    py = radius * disk_py;
    pz = radius * std::sqrt(U(1.0) - disk_px*disk_px - disk_py*disk_py);
}

template<typename T, typename U=typename detail::component<T>::type>
T hemisphere_cosine_concentric(U radius, U u1, U u2)
{
    T r;
    hemisphere_cosine_concentric(radius, u1, u2, r[0], r[1], r[2]);
    return r;
}

template<typename T, typename V, typename U=typename detail::component<T>::type>
T hemisphere_cosine_concentric(U radius, const V& u)
{
    T r;
    hemisphere_cosine_concentric(radius, u[0], u[1], r[0], r[1], r[2]);
    return r;
}

template<typename U>
void hemisphere_cosine_concentric(U radius, const U* u, U* p)
{
    hemisphere_cosine_concentric(radius, u[0], u[1], p[0], p[1], p[2]);
}

} // rsm
