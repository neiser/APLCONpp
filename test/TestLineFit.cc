#include "catch.hpp"

#include <APLCON.hpp>

#include "Value_t.h"

#include <algorithm>

using namespace std;

auto residuals = [] (const Value_t& a, const Value_t& b, const vector<Value_t>& x, const vector<Value_t>& y) {
    vector<double> residuals(y.size());
    transform(x.begin(), x.end(), y.begin(), residuals.begin(),
              [&a, &b] (const double& x_i, const double& y_i) {
        return a + b*x_i - y_i;
    });
    return residuals;
};

TEST_CASE("LineFit y-errors","") {
    vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
    vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

    for(auto& i : x)
        i.Fixed = true;

    // two unmeasured variables
    Value_t a{0,0};
    Value_t b{0,0};

    APLCON::Fitter<Value_t, Value_t, vector<Value_t>, vector<Value_t>> fitter;

    const auto& r = fitter.DoFit(a, b, x, y, residuals);

    CHECK(r.Status         == APLCON::Result_Status_t::Success);
    CHECK(r.ChiSquare      == Approx(67.44063286));
    CHECK(r.NDoF           == 2);
    CHECK(r.Probability    == Approx(0));
    CHECK(r.NIterations    == 2);
    CHECK(r.NFunctionCalls == 27);

    CHECK(a.Value() == Approx(0.0839608925));
    CHECK(a.Sigma() == Approx(0.0810352489));
    CHECK(b.Value() == Approx(0.9229444219));
    CHECK(b.Sigma() == Approx(0.0324371760));
}

TEST_CASE("LineFit xy-errors","") {
    vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
    vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

    // two unmeasured variables
    Value_t a{0,0};
    Value_t b{0,0};

    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = 10;

    APLCON::Fitter<Value_t, Value_t, vector<Value_t>, vector<Value_t>> fitter(settings);

    const auto& r = fitter.DoFit(a, b, x, y, residuals);

    CHECK(r.Status         == APLCON::Result_Status_t::Success);
    CHECK(r.ChiSquare      == Approx(18.4092764561));
    CHECK(r.NDoF           == 2);
    CHECK(r.Probability    == Approx(0.0001005718));
    CHECK(r.NIterations    == 8);
    CHECK(r.NFunctionCalls == 170);

    CHECK(a.Value() == Approx(-0.2338417299));
    CHECK(a.Sigma() == Approx( 0.2788696507));
    CHECK(b.Value() == Approx( 0.9786503645));
    CHECK(b.Sigma() == Approx( 0.1002995673));
}