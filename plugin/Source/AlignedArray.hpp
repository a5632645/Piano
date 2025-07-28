#pragma once
#include <cassert>
#include <corecrt_malloc.h>
#include <cstddef>
#include <memory>

template<class T>
class MySpan {
public:
    MySpan() = default;
    MySpan(T* ptr, size_t size) : ptr_(ptr), size_(size) {}
    T* Get() const { return ptr_; }
    T& operator[](size_t i) {
        assert(i < size_);
        return ptr_[i];
    }
private:
    T* ptr_{ nullptr };
    size_t size_{ 0 };
};

template<class T, size_t aligenment>
class AlignedArray {
public:
    ~AlignedArray() {
        Free();
    }

    void Free() {
        if (ptr_) {
            _aligned_free(ptr_);
            ptr_ = nullptr;
            size_ = 0;
        }
    }

    void Reset(size_t size) {
        Free();
        ptr_ = (T*)_aligned_malloc(size * sizeof(T), aligenment);
        std::fill_n(ptr_, size, T{});
        size_ = size;
    }

    T& operator[](size_t i) {
        assert(i < size_);
        return ptr_[i];
    }

    T operator[](size_t i) const {
        assert(i < size_);
        return ptr_[i];
    }

    T* Get() {
        return ptr_;
    }

    const T* Get() const {
        return ptr_;
    }

    MySpan<T> Span(size_t offset) const {
        assert(offset < size_);
        return {ptr_ + offset, size_ - offset};
    }

    MySpan<T> operator+(size_t offset) const {
        return Span(offset);
    }
private:
    T* ptr_{ nullptr };
    size_t size_{ 0 };
};
