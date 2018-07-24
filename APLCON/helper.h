#pragma once

#include <vector>

template<typename T>
class vec {
    std::vector<T> c;

public:

    void resize(size_t n) {
        c.resize(n);
    }

    void resize_and_reset(size_t n) {
        resize(n);
        std::fill(c.begin(), c.end(), T{});
    }

    T& operator[](size_t i) {
        return c[i-1];
    }

    const T& operator[](size_t i) const {
        return c[i-1];
    }

};

using vecd = vec<double>;
using veci = vec<int>;