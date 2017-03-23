#include "catch.hpp"

#include <APLCON.hpp>

#include "Value_t.h"

using namespace std;

TEST_CASE("Simple","") {

    APLCON::Fitter<Value_t, Value_t, Value_t, Value_t> fitter;

    Value_t BF_e_A  (0.105, 0.01);
    Value_t BF_e_B  (0.135, 0.03);
    Value_t BF_tau_A(0.095, 0.03);
    Value_t BF_tau_B( 0.14, 0.03);

    auto equal_firsttwo = [] (double a, double b, double, double) { return a - b; };
    auto equal_lasttwo = [] (double, double, double a, double b) { return a - b; };

    const auto& r = fitter.DoFit(BF_e_A, BF_e_B, BF_tau_A, BF_tau_B,
                                 equal_firsttwo, equal_lasttwo);

    CHECK(r.ChiSquare      == Approx(2.025));
    CHECK(r.NDoF           == 2);
    CHECK(r.Probability    == Approx(0.3633096218));

    CHECK(BF_e_A.Value()   == Approx(0.108));
    CHECK(BF_e_A.Sigma()   == Approx(0.009486833));
    CHECK(BF_e_A.Pull()    == Approx(0.9486832981));
    CHECK(BF_e_B.Value()   == Approx(0.108));
    CHECK(BF_e_B.Sigma()   == Approx(0.009486833));
    CHECK(BF_e_B.Pull()    == Approx(-0.9486832981));
    CHECK(BF_tau_A.Value() == Approx(0.1175));
    CHECK(BF_tau_A.Sigma() == Approx(0.0212132034));
    CHECK(BF_tau_A.Pull()  == Approx(1.0606601718));
    CHECK(BF_tau_B.Value() == Approx(0.1175));
    CHECK(BF_tau_B.Sigma() == Approx(0.0212132034));
    CHECK(BF_tau_B.Pull()  == Approx(-1.0606601718));

}
