#pragma once

#include <vector>
#include <ostream>

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

    friend std::ostream& operator<<(std::ostream& s, const vec& v) {
        s << '[';
        for(const auto& i : v.c) {
            s << i << ", ";
        }
        return s << "]";
    }
};

using vecd = vec<double>;
using veci = vec<int>;