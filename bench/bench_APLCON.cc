
#include "APLCON.hpp"

#include "benchmark/benchmark_api.h"

#include <ostream>
#include <vector>
#include <array>


struct X {

    constexpr X(double v, double s) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return N == APLCON::ValueIdx ? Value : Sigma;
    }

    friend std::ostream& operator<<(std::ostream& s, const X& o) {
        return s << "(" << o.Value << "," << o.Sigma << ")";
    }
};


static void BM_ErrorPropagation(benchmark::State& state) {
    auto a_and_b_is_c = [] (const X& a, const X& b, const X& c) {
        return c.Value - a.Value - b.Value;
    };
    while (state.KeepRunning()) {

        X a{10, 0.3};
        X b{20, 0.4};
        X c{ 0, 0.0};

        APLCON::Fitter<X, X, X> fitter;
        fitter.DoFit(a, b, c, a_and_b_is_c);

    }
}
BENCHMARK(BM_ErrorPropagation);

BENCHMARK_MAIN()

