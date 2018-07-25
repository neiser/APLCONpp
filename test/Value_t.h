#pragma once

#include "APLCON.hpp"
#include <ostream>

struct Value_t {

    constexpr Value_t(double v, double s) : V_S_P{v, s, 0}{}

    std::tuple<double,double,double> V_S_P; // short for Value, Sigma, Pull

    bool Fixed = false;
    bool Poisson = false;

    // fitter is able to access fields to fit them (last one pull is just written, actually)
    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return std::tie(std::get<N>(V_S_P));
    }

    // user can provide settings for each variabl
    template<size_t innerIdx>
    APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
        (void)outerIdx; // unused, provided to user method for completeness
        APLCON::Variable_Settings_t settings;
        if(Fixed)
            settings.StepSize = 0;
        if(Poisson)
            settings.Distribution = APLCON::Distribution_t::Poissonian;
        return settings;
    }

    // the following methods are just user convenience (not required by fitter)

    friend std::ostream& operator<<(std::ostream& s, const Value_t& o) {
        return s << "(" << o.Value() << "," << o.Sigma() << "," << o.Pull() << ")";
    }

    double Value() const {
        return std::get<APLCON::ValueIdx>(V_S_P);
    }

    double Sigma() const {
        return std::get<APLCON::SigmaIdx>(V_S_P);
    }

    double Pull() const {
        return std::get<APLCON::PullIdx>(V_S_P);
    }

    operator double() const {
        return std::get<APLCON::ValueIdx>(V_S_P);
    }
};