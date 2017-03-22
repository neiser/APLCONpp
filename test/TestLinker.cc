#include "catch.hpp"

#include <old/APLCON.hpp>

using namespace std;

TEST_CASE("Linker","") {

    struct BF_t {
        double A;
        double B;
    };

    auto linker = [] (BF_t& v) -> vector<double*> {
        return {addressof(v.A), addressof(v.B)}; // vector of two double*
    };

    auto equality_constraint = [] (vector<double> a) { return a[0] - a[1]; };

    {
        BF_t BF_e   = {0.105,0.135};
        BF_t BF_tau = {0.095,0.14};

        APLCON a("Linked variables");


        a.LinkVariable("BF_tau",linker(BF_tau), vector<double>{0.03,0.03});

        a.LinkVariable("BF_e", linker(BF_e), vector<double>{0.01,0.03});

        a.AddConstraint("BF_e_equal",   {"BF_e"},   equality_constraint);
        a.AddConstraint("BF_tau_equal", {"BF_tau"}, equality_constraint);

        a.DoFit();

        CHECK( BF_e.A   == Approx(0.108) );
        CHECK( BF_tau.A == Approx(0.1175) );
    }

    {
        APLCON b("Linked sigmas");

        BF_t BF_e   = {0.105, 0.135};
        BF_t BF_tau = {0.095, 0.14};
        BF_t BF_e_sigma = {0.01, 0.03};

        b.LinkVariable("BF_e",   linker(BF_e),   linker(BF_e_sigma));
        b.LinkVariable("BF_tau", linker(BF_tau), vector<double>{0.03,0.03});

        b.AddConstraint("BF_e_equal",   {"BF_e"},   equality_constraint);
        b.AddConstraint("BF_tau_equal", {"BF_tau"}, equality_constraint);

        b.DoFit();

        CHECK( BF_e.A       == Approx(0.108) );
        CHECK( BF_e_sigma.A == Approx(0.00949) );
    }
}
