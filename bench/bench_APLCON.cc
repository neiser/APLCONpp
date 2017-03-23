
#include "APLCON.hpp"

#include "benchmark/benchmark_api.h"

#include <ostream>
#include <vector>
#include <array>
#include <algorithm>

// borrow it from the tests...
#include "../test/Value_t.h"

using namespace std;

static void BM_ErrorPropagation(benchmark::State& state) {
    auto a_and_b_is_c = [] (const Value_t& a, const Value_t& b, const Value_t& c) {
        return c - a - b;
    };
    while (state.KeepRunning()) {

        {
            Value_t a{10, 0.3};
            Value_t b{20, 0.4};
            Value_t c{ 0, 0.0};

            APLCON::Fitter<Value_t, Value_t, Value_t> fitter;
            fitter.DoFit(a, b, c, a_and_b_is_c);
        }

        {
            Value_t a{10, 1};
            a.Poisson = true;
            Value_t b{20, 2};
            b.Poisson = true;
            Value_t c{ 0, 0.0};

            APLCON::Fitter<Value_t, Value_t, Value_t> fitter;
            fitter.DoFit(a, b, c, a_and_b_is_c);
        }
    }
}
BENCHMARK(BM_ErrorPropagation);

static void BM_LineFit(benchmark::State& state) {

    auto residuals = [] (const Value_t& a, const Value_t& b, const vector<Value_t>& x, const vector<Value_t>& y) {

        // we use again std::transform to calculate the residuals
        vector<double> residuals(y.size());
        transform(x.begin(), x.end(), y.begin(), residuals.begin(),
                  [&a, &b] (const double& x_i, const double& y_i) {
            return a + b*x_i - y_i;
        });
        return residuals;
    };

    while (state.KeepRunning()) {

        {
            APLCON::Fitter<Value_t, Value_t, vector<Value_t>, vector<Value_t>> fitter;

            vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
            vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

            for(auto& i : x)
                i.Fixed = true;

            // two unmeasured variables
            Value_t a{0,0};
            Value_t b{0,0};

            benchmark::DoNotOptimize(fitter.DoFit(a, b, x, y, residuals));

        }

        {
            APLCON::Fit_Settings_t settings;
            settings.MaxIterations = 10;

            APLCON::Fitter<Value_t, Value_t, vector<Value_t>, vector<Value_t>> fitter(settings);

            vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
            vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

            // two unmeasured variables
            Value_t a{0,0};
            Value_t b{0,0};

            benchmark::DoNotOptimize(fitter.DoFit(a, b, x, y, residuals));

        }
    }
}

BENCHMARK(BM_LineFit);


BENCHMARK_MAIN()

