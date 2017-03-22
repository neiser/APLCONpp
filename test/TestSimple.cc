#include "catch.hpp"

#include <old/APLCON.hpp>

using namespace std;

TEST_CASE("Simple","") {

    auto equality_constraint = [] (double a, double b) { return a - b; };

    APLCON a("WITHOUT lepton universality");

    a.AddMeasuredVariable("BF_e_A",   0.1050, 0.01);
    a.AddMeasuredVariable("BF_e_B",   0.135,  0.03);
    a.AddMeasuredVariable("BF_tau_A", 0.095,  0.03);
    a.AddMeasuredVariable("BF_tau_B", 0.14,   0.03);

    a.AddConstraint("BF_e_equal", {"BF_e_A", "BF_e_B"}, equality_constraint);
    a.AddConstraint("BF_tau_equal", {"BF_tau_A", "BF_tau_B"}, equality_constraint);

    const APLCON::Result_t& r = a.DoFit();

    CHECK(r.ChiSquare   == Approx(2.025));
    CHECK(r.NDoF        == 2);
    CHECK(r.Probability == Approx(0.3633096218));

    CHECK(r.Variables.at("BF_e_A").Value.After   == Approx(0.108));
    CHECK(r.Variables.at("BF_e_A").Sigma.After   == Approx(0.009486833));
    CHECK(r.Variables.at("BF_e_A").Pull          == Approx(0.9486832981));
    CHECK(r.Variables.at("BF_e_B").Value.After   == Approx(0.108));
    CHECK(r.Variables.at("BF_e_B").Sigma.After   == Approx(0.009486833));
    CHECK(r.Variables.at("BF_e_B").Pull          == Approx(-0.9486832981));
    CHECK(r.Variables.at("BF_tau_A").Value.After == Approx(0.1175));
    CHECK(r.Variables.at("BF_tau_A").Sigma.After == Approx(0.0212132034));
    CHECK(r.Variables.at("BF_tau_A").Pull        == Approx(1.0606601718));
    CHECK(r.Variables.at("BF_tau_B").Value.After == Approx(0.1175));
    CHECK(r.Variables.at("BF_tau_B").Sigma.After == Approx(0.0212132034));
    CHECK(r.Variables.at("BF_tau_B").Pull        == Approx(-1.0606601718));

    CHECK(r.Variables.at("BF_e_A").Covariances.After.at("BF_e_A")   == Approx(9E-05));
    CHECK(r.Variables.at("BF_e_A").Covariances.After.at("BF_e_B")   == Approx(9E-05));
    CHECK(r.Variables.at("BF_e_A").Covariances.After.at("BF_tau_A") == Approx(-0.0));
    CHECK(r.Variables.at("BF_e_A").Covariances.After.at("BF_tau_B") == Approx(-0.0));

    CHECK(r.Variables.at("BF_e_B").Covariances.After.at("BF_e_B")   == Approx(9E-05));
    CHECK(r.Variables.at("BF_e_B").Covariances.After.at("BF_tau_A") == Approx(-0.0));
    CHECK(r.Variables.at("BF_e_B").Covariances.After.at("BF_tau_B") == Approx(-0.0));

    CHECK(r.Variables.at("BF_tau_A").Covariances.After.at("BF_tau_A") == Approx(0.00045));
    CHECK(r.Variables.at("BF_tau_A").Covariances.After.at("BF_tau_B") == Approx(0.00045));
    CHECK(r.Variables.at("BF_tau_B").Covariances.After.at("BF_tau_B") == Approx(0.00045));
}
