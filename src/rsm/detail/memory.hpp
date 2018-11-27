/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdlib>

#ifndef RSM_DEFAULT_ALIGNMENT
#define RSM_DEFAULT_ALIGNMENT 16
#endif

namespace rsm {

struct allocator_t
{
    void* (*fn_alloc)(size_t, size_t, void*);
    void  (*fn_free)(void*, void*);
    void* user_data;
};

namespace detail {

inline const allocator_t& default_allocator()
{
    static allocator_t allocator = {
        [](size_t size, size_t alignment, void*) -> void* {
            // TODO: Handle alignment.
            (void)alignment;
            return std::malloc(size);
        },
        [](void* ptr, void*) {
            std::free(ptr);
        },
        nullptr,
    };
    return allocator;
}

template<typename T>
T* alloc(const allocator_t& allocator, size_t n)
{
    assert(allocator.fn_alloc);
    return reinterpret_cast<T*>(allocator.fn_alloc(sizeof(T) * n, RSM_DEFAULT_ALIGNMENT, allocator.user_data));
}

template<typename T>
void free(const allocator_t& allocator, T*& ptr)
{
    assert(allocator.fn_free);
    allocator.fn_free(ptr, allocator.user_data);
    ptr = nullptr;
}

} // detail
} // rsm
