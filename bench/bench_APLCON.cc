
#include "APLCON.hpp"

#include "benchmark/benchmark_api.h"

#include <ostream>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

using namespace std;

struct X {

    constexpr X(double v, double s) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return N == APLCON::ValueIdx ? std::tie(Value) : std::tie(Sigma);
    }

    friend std::ostream& operator<<(std::ostream& s, const X& o) {
        return s << "(" << o.Value << "," << o.Sigma << ")";
    }

    operator double() const {
        return Value;
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

static void BM_LineFit(benchmark::State& state) {

    auto residuals = [] (const X& a, const X& b, const vector<X>& x, const vector<X>& y) {

        // we use again std::transform to calculate the residuals
        vector<double> residuals(y.size());
        transform(x.begin(), x.end(), y.begin(), residuals.begin(),
                  [&a, &b] (const double& x_i, const double& y_i) {
            return a + b*x_i - y_i;
        });
        return residuals;
    };

    while (state.KeepRunning()) {

        vector<X> x{ {1,   0.2}, {2,   0.23}, {3,   0.16}, {4,   0.21} };
        vector<X> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

        // two unmeasured variables
        X a{0,0};
        X b{0,0};

        APLCON::Fitter<X, X, vector<X>, vector<X>> fitter;
        fitter.DoFit(a, b, x, y, residuals);
    }
}

BENCHMARK(BM_LineFit);


BENCHMARK_MAIN()

