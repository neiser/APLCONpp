#include "catch.hpp"

#include <old/APLCON.hpp>

using namespace std;

TEST_CASE("Very Simple", "") {

    APLCON a("Error propagation");

    a.AddMeasuredVariable("A", 10, 0.3);
    a.AddMeasuredVariable("B", 20, 0.4);

    a.AddUnmeasuredVariable("C");


    auto equality_constraint = [] (double a, double b, double c) { return c - a - b; };
    a.AddConstraint("A+B=C", {"A", "B", "C"}, equality_constraint);

    const APLCON::Result_t& ra = a.DoFit();

    REQUIRE(ra.Variables.at("C").Value.After == Approx(30.0));

    REQUIRE(ra.Variables.at("C").Sigma.After == Approx(0.5));

}
