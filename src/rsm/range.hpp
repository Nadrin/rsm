/*
 * rsm :: Random Sampling Mathematics
 * Copyright (c) 2018 Michał Siejak
 * Released under the MIT license; see LICENSE file for details.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <array>
#include <iterator>

namespace rsm {

template<unsigned int N, typename T>
class range_t
{
    static_assert(N > 0, "Number of dimensions must be greater than zero");
public:
    class iterator
    {
        friend class range_t;
    public:
        struct value_type {
            T* value;
            std::array<size_t, N> index;
        };
        struct value_type_wrapper {
            value_type_wrapper(value_type v) : value(v) {}
            value_type operator*() const { return value; }
            value_type value;
        };

        using reference = const value_type&;
        using pointer = const value_type*;
        using iterator_category = std::input_iterator_tag;
        using difference_type = void;

        reference operator*() const
        {
            return m_it;
        }
        pointer operator->() const
        {
            return &m_it;
        }
        friend bool operator==(const iterator& lhs, const iterator& rhs)
        {
            return lhs.m_it.value == rhs.m_it.value;
        }
        friend bool operator!=(const iterator& lhs, const iterator& rhs)
        {
            return lhs.m_it.value != rhs.m_it.value;
        }

        iterator& operator++()
        {
            for(unsigned int dim=0; dim<N; ++dim) {
                if(++m_it.index[dim] < m_size[dim]) {
                    m_it.value += m_stride;
                    return *this;
                }
                m_it.index[dim] = 0;
            }
            // Iterate past the buffer to produce the end iterator.
            m_it.value += m_stride;
            return *this;
        }
        value_type_wrapper operator++(int)
        {
            value_type_wrapper result(m_it);
            ++*this;
            return result;
        }

    private:
        iterator(const range_t* range, T* ptr)
            : m_size(range->m_size)
            , m_stride(range->m_stride)
        {
            m_it.value = ptr;
            m_it.index = {};
        }

        value_type m_it;
        std::array<size_t, N> m_size;
        uint16_t m_stride;
    };

    range_t(T* buffer, const std::array<size_t, N>& size, uint16_t stride=1)
        : m_buffer(buffer)
        , m_size(size)
        , m_stride(stride)
    {
        assert(m_stride > 0);
    }

    iterator begin() const
    {
        return iterator(this, m_buffer);
    }

    iterator end() const
    {
        return iterator(this, m_buffer + total_size() * m_stride);
    }

    size_t size(unsigned int dim) const
    {
        assert(dim < N);
        return m_size[dim];
    }

    size_t total_size() const
    {
        size_t s = m_size[0];
        for(size_t dim=1; dim<N; ++dim) {
            s *= m_size[dim];
        }
        return s;
    }

private:
    T* m_buffer;
    std::array<size_t, N> m_size;
    uint16_t m_stride;
};

template<unsigned int N, typename T>
range_t<N, T> range(T* buffer, const std::array<size_t, N>& size, uint16_t stride=1)
{
    return range_t<N, T>{buffer, size, stride};
}

} // rsm
