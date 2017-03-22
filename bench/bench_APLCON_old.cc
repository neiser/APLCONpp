#include "benchmark/benchmark_api.h"

#include <vector>
#include <iostream>

#include "APLCON.hpp"

using namespace std;

//static void BM_ErrorPropagation(benchmark::State& state) {

//    while(state.KeepRunning()) {

//        APLCON a("Error propagation");

//        a.AddMeasuredVariable("A", 10, 0.3);
//        a.AddMeasuredVariable("B", 20, 0.4);

//        a.AddUnmeasuredVariable("C");

//        auto equality_constraint = [] (double a, double b, double c) { return c - a - b; };
//        a.AddConstraint("A+B=C", {"A", "B", "C"}, equality_constraint);

//        // do the fit, obtain ra structure
//        const APLCON::Result_t& ra = a.DoFit();
//        cout << ra << endl;

//        APLCON b("Poissonian error propagation");

//        APLCON::Variable_Settings_t settings = APLCON::Variable_Settings_t::Default;
//        settings.Distribution = APLCON::Distribution_t::Poissonian;

//        b.AddMeasuredVariable("A", 10, 1, settings);
//        b.AddMeasuredVariable("B", 20, 2, settings);

//        b.AddUnmeasuredVariable("C");

//        b.AddConstraint("A+B=C", {"A", "B", "C"}, equality_constraint);

//        const APLCON::Result_t& rb = b.DoFit();
//        cout << rb << endl;

//    }
//}

//BENCHMARK(BM_ErrorPropagation);

static void BM_LineFit(benchmark::State& state) {

    struct data_t {
        vector<double> x;
        vector<double> sx; // errors in x
        vector<double> y;
        vector<double> sy; // errors in y
    };

    const data_t data{
        {   1,    2,    3,    4},
        { 0.2, 0.23, 0.16, 0.21},
        { 1.1, 1.95, 2.02, 3.98},
        {0.08, 0.04, 0.11, 0.07}
    };

    auto linker = [] (vector<double>& v) {
        vector<double*> v_p(v.size());
        transform(v.begin(), v.end(), v_p.begin(), addressof<double>);
        return v_p;
    };

    auto residuals = [] (const vector< vector<double> >& arg) {
        const double& a = arg[0][0];
        const double& b = arg[1][0];
        const vector<double>& x = arg[2];
        const vector<double>& y = arg[3];

        // we use again std::transform to calculate the residuals
        vector<double> residuals(y.size());
        transform(x.begin(), x.end(), y.begin(), residuals.begin(),
                  [&a, &b] (const double& x_i, const double& y_i) {
            return a + b*x_i - y_i;
        });
        return residuals;
    };

    while(state.KeepRunning()) {
        {
            APLCON f1("StraightLineFit");

            data_t data1 = data;

            APLCON::Variable_Settings_t var_settings = APLCON::Variable_Settings_t::Default;
            var_settings.StepSize = 0; // stepsize=0 means fixed variable
            f1.LinkVariable("x", linker(data1.x), vector<double>{0}, {var_settings}); // use x as fixed variables
            f1.LinkVariable("y", linker(data1.y), linker(data1.sy));

            f1.AddUnmeasuredVariable("a");
            f1.AddUnmeasuredVariable("b");

            f1.AddConstraint("residuals", vector<string>{"a", "b", "x", "y"}, residuals);

            benchmark::DoNotOptimize(f1.DoFit());
        }

        {

            APLCON f2("StraightLineFitWithXYErrors");

            data_t data2 = data;
            f2.LinkVariable("x", linker(data2.x), linker(data2.sx));
            f2.LinkVariable("y", linker(data2.y), linker(data2.sy));

            f2.AddUnmeasuredVariable("a");
            f2.AddUnmeasuredVariable("b");

            f2.AddConstraint("residuals", vector<string>{"a", "b", "x", "y"}, residuals);

            benchmark::DoNotOptimize(f2.DoFit());
        }
    }

}

BENCHMARK(BM_LineFit);


BENCHMARK_MAIN()
