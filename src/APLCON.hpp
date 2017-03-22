#pragma once

// link to the FORTRAN world
extern "C" {
#include "wrapper/APLCON.h"
}

// some template in extra header file
#include "detail/APLCON_hpp.hpp"
#include <cmath>  // for sqrt
#include <limits>
#include <functional>

namespace APLCON {

// use those constants in your classes to
// identify what is request in the call to "linkFitter"
static constexpr auto ValueIdx = 0;
static constexpr auto SigmaIdx = 1;
static constexpr auto PullIdx  = 2;

static constexpr auto Inf = std::numeric_limits<double>::infinity();
static constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();

/**
 * @brief The Distribution_t enum
 * @see Variable_Settings_t
 */
enum class Distribution_t {
    Gaussian, /**< Gaussian distributed variable (default) */
    Poissonian, /**< Poissonian distributed variable */
    LogNormal, /**< Ratios are lognormal distributed */
    SquareRoot /**< SquareRoot transformation */
};

/**
 * @brief The Limit_t struct defines upper and lower limits
 */
struct Limit_t {
    double Low;
    double High;
};

/**
 * @brief The Variable_Settings_t struct contains settings per variable
 */
struct Variable_Settings_t {
    Distribution_t Distribution = Distribution_t::Gaussian;
    Limit_t Limit{-Inf, +Inf};
    double StepSize{NaN};
};

/**
 * @brief The Result_Status_t enum encodes APLCON's status after fit
 * @note the order must correspond to APLCON's status number (see APLCON README)
 */
enum class Result_Status_t : int {
    Success, /**< Fit was successful, result is meaningful */
    NoConvergence, /**< No convergence reached */
    TooManyIterations, /**< Too many iterations needed with no convergence */
    UnphysicalValues, /**< Unphysical values encountered during fit */
    NegativeDoF, /**< Negative degrees of freedom */
    OutOfMemory, /**< Not sufficient memory for fit */
    _Unknown // default in Result_t, also used to count items
};

/**
 * @brief The Result_t struct contains after the fit all information about it.
 */
struct Result_t {
    Result_Status_t Status{Result_Status_t::_Unknown};
    double ChiSquare{NaN};
    int NDoF{-1};
    double Probability{NaN};
    int NIterations{-1};
    int NFunctionCalls{-1};
};

/**
* @brief The Fit_Settings_t struct. See APLCON itself for details.
*/
struct Fit_Settings_t {
    int DebugLevel{0}; // no debug output
    int MaxIterations{-1}; // default
    double ConstraintAccuracy{NaN}; //
    double Chi2Accuracy{NaN};
    double MeasuredStepSizeFactor{NaN};
    double UnmeasuredStepSizeFactor{NaN};
    double MinimalStepSizeFactor{NaN};
};

/**
 * @brief The Fitter carries out the fit with DoFit
 */
template<class ... Vars>
struct Fitter {

    explicit Fitter(const Fit_Settings_t& settings = {}) :
        FitSettings(settings) {}

    template<class ... Constraints>
    Result_t DoFit(Vars&... vars, Constraints&&... constraints) {

        using namespace detail;

        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nv>());
        {
            constexpr auto innerDim_sigmas = getInnerDim<SigmaIdx>(build_indices<Nv>());
            static_assert(compare_array(innerDim, innerDim_sigmas), "Method linkFitter returns unequal number of values and sigmas.");
        }

        OuterDim = getOuterDim(vars...);

        const auto nVar = sum_of_array(prod_of_array(OuterDim, innerDim));

        // fill the initial values
        {
            X.resize(nVar);
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (double& v, const double& t) { v = t; }, vars...);
        }

        // fill the initial uncertainties / sigmas
        {
            V.resize((nVar*nVar+nVar)/2);
            linker_diagonal_t linker(V.begin());
            callLinkFitter<SigmaIdx>(linker, [] (double& v, const double& t) { v = t*t; }, vars...);
        }


        // dispatch via AllocF to ensure optimal constraint storage handling
        // (and then reach RunAPLCON)
        constexpr auto isConstexpr = sum_of(isConstraintsSizeConstexpr<Constraints>()...)==0;
        return AllocF(std::enable_if<isConstexpr>(), vars..., std::forward<Constraints>(constraints)...);
    }

private:
    Fit_Settings_t FitSettings;

    static constexpr auto Nv = sizeof...(Vars);
    using idx_array_t = std::array<std::size_t, Nv>;

    // keep storage here to avoid re-allocations
    std::vector<double> X;
    std::vector<double> V;
    std::vector<double> F_dynamic; // only used if non-constexpr constraints size

    idx_array_t OuterDim;

    template<class ... Constraints>
    Result_t AllocF(std::enable_if<true>, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;
        constexpr auto nConstraints = sum_of(constraint_test<Constraints(Vars...)>().getN()...);
        // alloc on stack
        std::array<double, nConstraints> F_static;
        return RunAPLCON(F_static, vars..., std::forward<Constraints>(constraints)...);
    }

    template<class ... Constraints>
    Result_t AllocF(std::enable_if<false>, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;
        const auto nConstraints = sum_of(getConstraintDim<Constraints>(std::forward<Constraints>(constraints), vars...)...);
        // alloc on heap
        F_dynamic.resize(nConstraints);
        return RunAPLCON(F_dynamic, vars..., std::forward<Constraints>(constraints)...);
    }

    template<class F_t, class ... Constraints>
    Result_t RunAPLCON(F_t& F, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;

        c_aplcon_aplcon(X.size(), F.size());
        c_aplcon_aprint(6, FitSettings.DebugLevel); // default output on LUNP 6 (STDOUT)

        call_if_set(c_aplcon_apdeps,   FitSettings.ConstraintAccuracy);
        call_if_set(c_aplcon_apepschi, FitSettings.Chi2Accuracy);
        call_if_set(c_aplcon_apiter,   FitSettings.MaxIterations);
        call_if_set(c_aplcon_apderf,   FitSettings.MeasuredStepSizeFactor);
        call_if_set(c_aplcon_apderu,   FitSettings.UnmeasuredStepSizeFactor);
        call_if_set(c_aplcon_apdlow,   FitSettings.MinimalStepSizeFactor);

        // still don't use simple loop,
        // as we can check at compile time if something to do
        setVarSettings(ctie(vars...));

        // the main convergence loop
        int aplcon_ret = -1;
        do {

            // compute F from constraints
            callConstraints(F.begin(), vars..., std::forward<Constraints>(constraints)...);

            // call APLCON iteration
            c_aplcon_aploop(X.data(), V.data(), F.data(), &aplcon_ret);

            // copy X into vars, for constraint evaluation
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (const double& v, double& t) { t = v; }, vars...);

        }
        while(aplcon_ret<0);

        {
            // copy X into vars, one last time
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (const double& v, double& t) { t = v; }, vars...);
        }

        {
            // copy back the sigmas
            // don't forget that sigmas are the sqrt's of the covariance diagonal
            linker_diagonal_t linker(V.begin());
            callLinkFitter<SigmaIdx>(linker, [] (const double& v, double& t) { t = std::sqrt(v); }, vars...);
        }

        {
            // get the pulls from APLCON
            std::vector<double> pulls(X.size());
            c_aplcon_appull(pulls.data());

            // copy pulls into vars
            linker_linear_t linker(pulls.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<PullIdx>(linker, [] (const double& v, double& t) { t = v; }, vars...);
        }

        // now retrieve "everything" from APLCON
        Result_t result;
        result.Status = static_cast<Result_Status_t>(aplcon_ret);

        // retrieve some info about the fit (directly copy to struct field if possible)
        // chndpv and apstat both return the resulting chi2,
        // but the latter returns it with double precision
        float chi2, pval;
        c_aplcon_chndpv(&chi2,&result.NDoF,&pval);
        result.Probability = pval;
        c_aplcon_apstat(&result.ChiSquare, &result.NFunctionCalls, &result.NIterations);

        return result;
    }

    /// setVarSettings

    template<std::size_t I = 0, class... Vars_>
    typename std::enable_if<I < sizeof...(Vars), void>::type
    setVarSettings(const std::tuple<Vars...>& t) noexcept {
        using namespace detail;
        using Var = typename std::tuple_element<I, std::tuple<Vars...>>::type;
        using uVar = typename decay_stl_cont<Var>::type; // look inside stl container types for its contained element type (=value_type)
        setVarSettings<Var, I>(std::enable_if<has_getFitterSettings<uVar, Variable_Settings_t>::value>(), // dispatch has method
                               std::get<I>(t));
        setVarSettings<I+1, Vars...>(t);
    }

    // index counting recursion end
    template<std::size_t I = 0, class... Vars_>
    typename std::enable_if<I == sizeof...(Vars), void>::type
    setVarSettings(const std::tuple<Vars...>&) noexcept {}

    template<class Var, size_t I>
    void
    setVarSettings(std::enable_if<true>, const Var& var) noexcept {
        using namespace detail;
        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nv>());
        constexpr auto nInner = std::get<I>(innerDim);

        // nOuter is one if Var is not stl_container
        const auto nOuter = std::get<I>(OuterDim);

        // build_indices runs from 0...I-1, fits perfectly to sum up the offset
        const auto offset = sum_of_array(prod_of_array_impl(innerDim, OuterDim, build_indices<I>()));


        for(auto iOuter=0u;iOuter<nOuter;iOuter++) {
            const auto& var_settings = getVarSettings(std::enable_if<is_stl_cont<Var>::value>(),
                                                      build_indices<nInner>(), var, iOuter);

            for(auto iInner=0;iInner<nInner;iInner++) {
                // remember that FORTRAN/APLCON starts counting at 1
                const auto varidx = offset + iOuter*nInner + iInner + 1;
                const auto& s = var_settings[iInner];

                using namespace std::placeholders;

                // setup APLCON variable specific things
                switch (s.Distribution) {
                case APLCON::Distribution_t::Gaussian:
                    // thats the APLCON default, nothing must be called
                    break;
                case APLCON::Distribution_t::Poissonian:
                    c_aplcon_apoiss(varidx);
                    break;
                case APLCON::Distribution_t::LogNormal:
                    c_aplcon_aplogn(varidx);
                    break;
                case APLCON::Distribution_t::SquareRoot:
                    c_aplcon_apsqrt(varidx);
                    break;
                    // APLCON exposes even more transformations (see wrapper),
                    // but they're not mentioned in the README...
                default:
                    break;
                }

                call_if_set(std::bind(c_aplcon_aplimt, varidx, _1, _2), s.Limit.Low, s.Limit.High);
                call_if_set(std::bind(c_aplcon_apstep, varidx, _1), s.StepSize);
            }
        }
    }

    template<size_t... Idx, class Var>
    static constexpr std::array<Variable_Settings_t, sizeof...(Idx)>
    getVarSettings(std::enable_if<true>, detail::indices<Idx...>, const Var& var, size_t iOuter) {
        // Var is STL container, so ask contained item for fitter settings
        return {  std::next(var.begin(), iOuter)->template getFitterSettings<Idx>(iOuter)...  };
    }

    template<size_t... Idx, class Var>
    static constexpr std::array<Variable_Settings_t, sizeof...(Idx)>
    getVarSettings(std::enable_if<false>, detail::indices<Idx...>, const Var& var, size_t iOuter) {
        // Var is not a container, just call getFitterSettings then
        return {  var.template getFitterSettings<Idx>(iOuter)...  };
    }

    // do nothing if no getFitterSettings implemented
    template<class Var, size_t I>
    static void
    setVarSettings(std::enable_if<false>, const Var&) noexcept {}

    /// callConstraints

    template<class It, class Constraint, class... Constraints>
    static void callConstraints(It it, const Vars&... vars, Constraint&& constraint, Constraints&&... constraints) noexcept {
        callConstraint<It>(it, std::forward<Constraint>(constraint), vars...);
        callConstraints(it, vars..., std::forward<Constraints>(constraints)...);
    }

    template<class It>
    static void callConstraints(It, const Vars&...)  noexcept {}

    template<class It, class Constraint>
    static typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_multi>::value,
    void>::type
    callConstraint(It& it, Constraint&& constraint, const Vars&... vars) noexcept {
        // Making F_ constexpr is probably too limiting for constraints (think of TreeFitter with rather complex lambdas)
        const auto& F_ = std::forward<Constraint>(constraint)(vars...);
        it = std::copy(F_.begin(), F_.end(), it);
    }

    template<class It, class Constraint>
    static typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_single>::value,
    void>::type
    callConstraint(It& it, Constraint&& constraint, const Vars&... vars) noexcept {
        // single scalar is just copied to current position
        *it++ = std::forward<Constraint>(constraint)(vars...);
    }

    /// callLinkFitter

    template<size_t N, class Linker, class Func, class Var, class... Vars_>
    static void callLinkFitter(Linker&& linker, Func&& f, Var&& var, Vars_&&... vars) noexcept {
        callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f),
                          std::enable_if<detail::is_stl_cont<Var>::value>(), std::forward<Var>(var));
        callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f), std::forward<Vars_>(vars)...);
    }

    template<size_t N, class Linker, class Func>
    static void callLinkFitter(Linker&&, Func&&) noexcept {}

    template<size_t N, class Linker, class Func, class Var>
    static void callLinkFitter(Linker&& linker, Func&& f, std::enable_if<false>, Var&& var) noexcept {
        std::forward<Linker>(linker)(std::forward<Func>(f), std::forward<Var>(var).template linkFitter<N>());
    }

    template<size_t N, class Linker, class Func, class Var>
    static void callLinkFitter(Linker&& linker, Func&& f, std::enable_if<true>, Var&& container) noexcept {
        for(auto& var : std::forward<Var>(container))
            callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f), std::enable_if<false>(), var);
    }

    /// getOuterDim

    template<class... Vars_>
    static idx_array_t
    getOuterDim(Vars_&&... vars) noexcept {
        using namespace detail;
        // expand and dispatch
        return {getOuterDim<Vars_>(
                        std::enable_if<is_stl_cont<Vars>::value>(), // dispatch
                        std::forward<Vars_>(vars)
                        )...}; // expand Vars_
    }

    template<class Var>
    static size_t
    getOuterDim(std::enable_if<true>, Var&& var) noexcept {
        // no other way than calling non-constexpr .size()
        return std::forward<Var>(var).size();
    }

    template<class Var>
    static constexpr size_t
    getOuterDim(std::enable_if<false>, Var&&) noexcept {
        return 1;
    }

    /// getInnerDim

    template<size_t N, size_t... Idx>
    static constexpr idx_array_t
    getInnerDim(detail::indices<Idx...>) noexcept {
        // dispatch depending if Var in parameter pack is container or not
        // if yes, inspect value_type. if no, inspect inner Var itself
        return {getInnerDim<N, Vars>(std::enable_if<detail::is_stl_cont<Vars>::value>())...};
    }

    template<size_t N, class Var>
    static constexpr std::size_t
    getInnerDim(std::enable_if<false>) noexcept {
        // inspect the templated function
        using f_t = decltype(&Var::template linkFitter<N>);
        using r_t = typename std::result_of<f_t(Var)>::type;
        return std::tuple_size<r_t>::value;

    }

    template<size_t N, class Var>
    static constexpr size_t
    getInnerDim(std::enable_if<true>) noexcept {
        return getInnerDim<N, typename Var::value_type>(std::enable_if<false>());
    }


    /// getConstraintDim

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_constsize>::value,
    std::size_t>::type
    getConstraintDim(Constraint&&, const Vars&...) noexcept {
        return detail::constraint_test<Constraint(Vars...)>().getN();
    }

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_container>::value,
    std::size_t>::type
    getConstraintDim(Constraint&& constraint, const Vars&... vars) noexcept {
        // call the constraint with the given vars (no other way to figure it out)
        return std::forward<Constraint>(constraint)(vars...).size();
    }

    /// isConstraintsSizeConstexpr

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_constsize>::value,
    std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 0;
    }

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    detail::constraint_test<Constraint(Vars...), detail::c_is_container>::value,
    std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 1;
    }
};

} // namespace APLCON