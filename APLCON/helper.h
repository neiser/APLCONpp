#pragma once

#include <vector>

template<typename T>
class vec {
    std::vector<T> c;

public:

    void resize(size_t n) {
        c.resize(n);
    }

    T& operator[](size_t i) {
        return c[i-1];
    }
};

using vecd = vec<double>;