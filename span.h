/*
 * span.h
 *
 *  Created on: Jun 18, 2025
 *      Author: sergey
 *
 *  Простой спан для вектора под С++03 (без проверок на допустимые типы)
 *
 */

#ifndef SPAN_H_
#define SPAN_H_

#if __cplusplus >= 202002L
#include <span>
#else

#include <cstddef>
#include <vector>
#include <type_traits>

namespace std {

namespace detail {

template <typename E>
struct span_storage {
    span_storage()
        : ptr(NULL), size(0) {}

    span_storage(E* p_ptr, size_t p_size)
        : ptr(p_ptr), size(p_size) {}

    E* ptr;
    size_t size;
};

}  // namespace detail

template <typename ElementType>
class span
{
    typedef detail::span_storage<ElementType> storage_type;

public:
    typedef ElementType element_type;
    typedef typename std::remove_cv<ElementType>::type value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef element_type* pointer;
    typedef const element_type* const_pointer;
    typedef element_type& reference;
    typedef const element_type& const_reference;
    typedef pointer iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    span()
    {}

    template <typename _It>
    span(_It first, size_type count)
        : storage_(&(*first), count)
    {}

    template <typename _It>
    span(_It first, _It last)
        : storage_(&(*first), std::distance(first, last))
    {}

    template <typename OtherElementType>
    span(const std::vector<OtherElementType>& cont)
        : storage_(&cont[0], cont.size())
    {}

    template <typename OtherElementType>
    span(std::vector<OtherElementType>& cont)
        : storage_(&cont[0], cont.size())
    {}

    template <typename OtherElementType>
    span(const span<OtherElementType>& other)
         : storage_(other.data(), other.size())
    {}

    span<element_type> first(size_type count) const
    {
        return span<element_type>(data(), count);
    }

    span<element_type> last(size_type count) const
    {
        return span<element_type>(data() + (size() - count), count);
    }

    size_type size() const
    {
        return storage_.size;
    }

    size_type size_bytes() const
    {
        return size() * sizeof(element_type);
    }

    bool empty() const
    {
        return size() == 0;
    }

    reference operator[](size_type idx) const
    {
        return *(data() + idx);
    }

    reference front() const
    {
        return *data();
    }

    reference back() const
    {
        return *(data() + (size() - 1));
    }

    pointer data() const
    {
        return storage_.ptr;
    }

    iterator begin() const
    {
        return data();
    }

    iterator end() const
    {
        return data() + size();
    }

    reverse_iterator rbegin() const
    {
        return reverse_iterator(end());
    }

    reverse_iterator rend() const
    {
        return reverse_iterator(begin());
    }

private:
    storage_type storage_;
};

} // namespace std {
#endif /* __cplusplus >= 202002L */

#endif /* SPAN_H_ */
