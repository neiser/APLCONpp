#include "catch.hpp"

#include <APLCON.hpp>

#include "Value_t.h"

using namespace std;

void dofit1()
{
    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-2; // needed to get it working

    Value_t N_fit(1328.90,104.124);
    Value_t N_mcreco(101430.,321.205);
    Value_t N_mcgen(263085.,522.916);
    Value_t N_effcorr(0,0);
    APLCON::Fitter<Value_t, Value_t, Value_t, Value_t> fitter(fit_settings);
    auto r = fitter.DoFit(N_fit, N_mcreco, N_mcgen, N_effcorr,
                                 [] (const Value_t& N_fit, const Value_t& N_mcreco,
                                 const Value_t& N_mcgen, const Value_t& N_effcorr) {
        return N_effcorr.Value() - N_fit.Value() * N_mcgen.Value() / N_mcreco.Value();
    });

    REQUIRE(r.Status == APLCON::Result_Status_t::Success);
    REQUIRE(N_effcorr.Sigma() == Approx(270.379889145));
}

void dofit2()
{
    APLCON::Fit_Settings_t fit_settings;
    fit_settings.ConstraintAccuracy = 1e-2; // needed to get it working

    std::vector<Value_t> Ns_fit{
        {1376.65, 0059.7718},
        {1328.90, 0104.124},
    };

    Value_t Nsum(0,0); // sigma=0 means unmeasured

    APLCON::Fitter<std::vector<Value_t>, Value_t> fitter(fit_settings);
    unsigned call = 0;
    auto r = fitter.DoFit(Ns_fit, Nsum, [&call] (
                          const vector<Value_t>& N, const Value_t& Nsum) {
        call++;
        double sum = 0.0;
        for(auto& n : N) {
            sum += n.Value();
        }
        return sum - Nsum.Value();
    });

    REQUIRE(call == 15);
    REQUIRE(r.Status == APLCON::Result_Status_t::Success);
    REQUIRE(Nsum.Sigma() == Approx(120.0602992302));
}

TEST_CASE("Uninitialized mem access","") {
    // although this test might seem trivial,
    // it helped debugging a nasty access of uninitialized memory acess
    // in the ENTRY call of DUMINV in "condutil.F".
    // We used this code here with a lot of WRITE statements inside APLCON
    // and compared the output with the help of `meld`
    // before fix, dofit2() only worked if dofit1() was NOT called
    dofit2();
    dofit1();
    dofit2();
}

