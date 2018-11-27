/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cstdint>

namespace rsm {

using options_t = uint8_t;

namespace opt {

enum option_bits {
    none    = 0,
    jitter  = 1,
    shuffle = 2,
};

} // opt
} // rsm
