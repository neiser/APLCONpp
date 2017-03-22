
#include "APLCON.hpp"

#include <ostream>
#include <iostream>
using namespace std;

struct X {

    constexpr X(double v, double s) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
//        return N == APLCON::ValueIdx ? std::tie(Value,Value) : std::tie(Sigma,Sigma);
        return N == APLCON::ValueIdx ? std::tie(Value) : std::tie(Sigma);
    }

    template<size_t N>
    APLCON::Variable_Settings_t getFitterSettings() const noexcept {
        APLCON::Variable_Settings_t settings;

        return settings;
    }

    friend std::ostream& operator<<(std::ostream& s, const X& o) {
        return s << "(" << o.Value << "," << o.Sigma << ")";
    }

    operator double() const {
        return Value;
    }
};

int main() {
    auto a_and_b_is_c = [] (const X& a, const X& b, const X& c) {
        return c.Value - a.Value - b.Value;
    };

    X a{10, 0.3};
    X b{20, 0.4};
    X c{ 0, 0.0};

    APLCON::Fitter<X, X, X> fitter;
    fitter.DoFit(a, b, c, a_and_b_is_c);

    cout << c << endl;

}
