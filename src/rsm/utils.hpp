/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstdint>
#include <utility>

#include "next.hpp"

namespace rsm {

template<typename Generator, typename T>
void shuffle(Generator& generator, T* begin, T* end)
{
    assert(end >= begin);
    const uint32_t count = end - begin;
    for(uint32_t i=0; i<count; ++i) {
        uint32_t r = next(generator, i, count);
        std::swap(begin[i], begin[r]);
    }
}

template<size_t BlockSize, typename Generator, typename T>
void shuffle(Generator& generator, T* begin, T* end)
{
    static_assert(BlockSize > 0, "BlockSize must not be zero");
    assert(end >= begin);
    const uint32_t count = (end - begin) / BlockSize;
    for(uint32_t i=0; i<count; ++i) {
        uint32_t r = next(generator, i, count);
        for(size_t j=0; j<BlockSize; ++j) {
            std::swap(begin[BlockSize*i+j], begin[BlockSize*r+j]);
        }
    }
}

template<size_t BlockSize, typename Generator, typename T>
void shuffle_inner(Generator& generator, T* begin, T* end)
{
    static_assert(BlockSize > 0, "BlockSize must not be zero");
    assert(end >= begin);
    const uint32_t count = (end - begin) / BlockSize;
    for(size_t j=0; j<BlockSize; ++j) {
        for(uint32_t i=0; i<count; ++i) {
            uint32_t r = next(generator, i, count);
            std::swap(begin[BlockSize*i+j], begin[BlockSize*r+j]);
        }
    }
}

} // rsm
