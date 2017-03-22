
#include "APLCON.hpp"

#include <ostream>
#include <iostream>
#include <algorithm>

using namespace std;

struct Value_t {

    constexpr Value_t(double v, double s, bool fixed = false) : Value(v), Sigma(s), Fixed(fixed) {}

    double Value;
    double Sigma;
    bool   Fixed;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return N == APLCON::ValueIdx ? std::tie(Value) : std::tie(Sigma);
//        return N == APLCON::ValueIdx ? std::tie(Value) : std::tie(Sigma);
    }

    template<size_t innerIdx>
    APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
        (void)outerIdx; // unused, provided to user struct for completeness
        APLCON::Variable_Settings_t settings;
        if(Fixed)
            settings.StepSize = 0;
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
    auto residuals = [] (const Value_t& a, const Value_t& b, const vector<Value_t>& x, const vector<Value_t>& y) {

        // we use again std::transform to calculate the residuals
        vector<double> residuals(y.size());
        transform(x.begin(), x.end(), y.begin(), residuals.begin(),
                  [&a, &b] (const double& x_i, const double& y_i) {
            return a + b*x_i - y_i;
        });
        return residuals;
    };

    APLCON::Fit_Settings_t settings;
    settings.DebugLevel = 5;
    settings.MaxIterations = 10;

    APLCON::Fitter<Value_t, Value_t, vector<Value_t>, vector<Value_t>> fitter{settings};

    {

        vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
        vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

        for(auto& i : x)
            i.Fixed = true;

        // two unmeasured variables
        Value_t a{0,0};
        Value_t b{0,0};

        fitter.DoFit(a, b, x, y, residuals);

        cout << "y-errors only: a=" << a << " b=" << b << endl;
    }

    cout << endl << endl << "<<<<<<<<<<<<<<<<<<<<<<<< " << endl << endl << endl;

    {
        vector<Value_t> x{ {1,  0.2},   {2,  0.23}, {3,   0.16}, {4,   0.21} };
        vector<Value_t> y{ {1.1,0.08}, {1.95,0.04}, {2.02,0.11}, {3.98,0.07} };

        for(auto& i : x)
            i.Fixed = false;

        // two unmeasured variables
        Value_t a{0,0};
        Value_t b{0,0};

        fitter.DoFit(a, b, x, y, residuals);

        cout << "xy-errors:     a=" << a << " b=" << b << endl;

    }
}
