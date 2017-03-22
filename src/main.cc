
#include "APLCON.hpp"

#include <ostream>
#include <iostream>
using namespace std;

struct Value_t {

    constexpr Value_t(double v, double s) : Value(v), Sigma(s) {}

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

    friend std::ostream& operator<<(std::ostream& s, const Value_t& o) {
        return s << "(" << o.Value << "," << o.Sigma << ")";
    }

    operator double() const {
        return Value;
    }
};

int main() {
    auto a_and_b_is_c = [] (const Value_t& a, const Value_t& b, const Value_t& c) {
        return c.Value - a.Value - b.Value;
    };

    Value_t a{10, 0.3};
    Value_t b{20, 0.4};
    Value_t c{ 0, 0.0};

    APLCON::Fitter<Value_t, Value_t, Value_t> fitter;
    fitter.DoFit(a, b, c, a_and_b_is_c);

    cout << c << endl;

}
