#include "catch.hpp"

#include <APLCON.hpp>

#include "Value_t.h"

using namespace std;

double a_and_b_is_c(const Value_t& a, const Value_t& b, const Value_t& c) {
    return c - a - b;
}

TEST_CASE("Very Simple 1", "") {

    APLCON::Fitter<Value_t, Value_t, Value_t> fitter;

    Value_t a{10, 0.3};
    Value_t b{20, 0.4};
    Value_t c{ 0, 0.0};

    const auto& result =  fitter.DoFit(a, b, c, a_and_b_is_c);

    REQUIRE(result.Status == APLCON::Result_Status_t::Success);

    REQUIRE(std::get<APLCON::ValueIdx>(c.V_S_P) == Approx(30.0));
    REQUIRE(std::get<APLCON::SigmaIdx>(c.V_S_P) == Approx(0.5));
}

TEST_CASE("Very Simple 2", "") {

    APLCON::Fitter<Value_t, Value_t, Value_t> fitter;

    // uncertainties don't matter obviously for Poisson
    Value_t a{10, 0};
    a.Poisson = true;
    Value_t b{20, 0};
    b.Poisson = true;
    Value_t c{ 0, 0.0};

    const auto& result =  fitter.DoFit(a, b, c, a_and_b_is_c);

    REQUIRE(result.Status == APLCON::Result_Status_t::Success);

    REQUIRE(std::get<APLCON::ValueIdx>(c.V_S_P) == Approx(30.0));

    REQUIRE(std::get<APLCON::SigmaIdx>(c.V_S_P) == Approx(sqrt(30.0)).epsilon(0.01));

}
