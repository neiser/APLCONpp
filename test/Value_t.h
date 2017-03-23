#pragma once

struct Value_t {

    constexpr Value_t(double v, double s) : V_S_P{v, s, 0}{}

    std::tuple<double,double,double> V_S_P;

    bool Fixed = false;
    bool Poisson = false;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        return std::tie(std::get<N>(V_S_P));
    }

    template<size_t innerIdx>
    APLCON::Variable_Settings_t getFitterSettings(size_t outerIdx) const noexcept {
        (void)outerIdx; // unused, provided to user struct for completeness
        APLCON::Variable_Settings_t settings;
        if(Fixed)
            settings.StepSize = 0;
        if(Poisson)
            settings.Distribution = APLCON::Distribution_t::Poissonian;
        return settings;
    }

    friend std::ostream& operator<<(std::ostream& s, const Value_t& o) {
        return s << "(" << std::get<APLCON::ValueIdx>(o.V_S_P)
                 << "," << std::get<APLCON::SigmaIdx>(o.V_S_P)
                 << "," << std::get<APLCON::PullIdx>(o.V_S_P)
                 << ")";
    }

    operator double() const {
        return std::get<APLCON::ValueIdx>(V_S_P);
    }
};