/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <cmath>

namespace rsm {

template<typename T>
constexpr T balance_heuristic(size_t n_f, T f_pdf, size_t n_g, T g_pdf)
{
    return (n_f * f_pdf) / (n_f * f_pdf + n_g * g_pdf);
}

template<typename T, int B=2>
constexpr T power_heuristic(size_t n_f, T f_pdf, size_t n_g, T g_pdf)
{
    T pow_f = std::pow(n_f * f_pdf, B);
    T pow_g = std::pow(n_g * g_pdf, B);
    return pow_f / (pow_f + pow_g);
}

} // rsm
