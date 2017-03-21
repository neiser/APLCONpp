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

// use those constants in your classes to
// identify what is request in the call to "linkFitter"

namespace APLCON {

static constexpr auto ValueIdx = 0;
static constexpr auto SigmaIdx = 1;

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
 * @brief The Fitter carries out the fit with DoFit
 */
template<class ... Vars>
struct Fitter {

    template<class ... Constraints>
    void DoFit(Vars&... vars, Constraints&&... constraints) {

        using namespace detail;

        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nv>());
        {
            constexpr auto innerDim_sigmas = getInnerDim<SigmaIdx>(build_indices<Nv>());
            static_assert(compare_array(innerDim, innerDim_sigmas), "The implementation of link_fitter returns unequal number of values and sigmas.");
        }
        const auto nVar = getNvars(innerDim, build_indices<Nv>(), vars...);

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
        AllocF(std::enable_if<isConstexpr>(), vars..., std::forward<Constraints>(constraints)...);

    }

private:

    // keep storage here to avoid re-allocations
    std::vector<double> X;
    std::vector<double> V;
    std::vector<double> F_dynamic; // only used if non-constexpr constraints size

    static constexpr auto Nv = sizeof...(Vars);
    using idx_array_t = std::array<std::size_t, Nv>;

    template<class ... Constraints>
    void AllocF(std::enable_if<true>, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;
        constexpr auto nConstraints = sum_of(constraint_test<Constraints(Vars...)>().getN()...);
        // alloc on stack
        std::array<double, nConstraints> F_static;
        RunAPLCON(F_static, vars..., std::forward<Constraints>(constraints)...);
    }

    template<class ... Constraints>
    void AllocF(std::enable_if<false>, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;
        const auto nConstraints = sum_of(getConstraintDim<Constraints>(std::forward<Constraints>(constraints), vars...)...);
        // alloc on heap
        F_dynamic.resize(nConstraints);
        RunAPLCON(F_dynamic, vars..., std::forward<Constraints>(constraints)...);
    }

    template<class F_t, class ... Constraints>
    void RunAPLCON(F_t& F, Vars&... vars, Constraints&&... constraints) {
        using namespace detail;

        c_aplcon_aplcon(X.size(), F.size());

        c_aplcon_aprint(6, 0); // default output on LUNP 6 (STDOUT)

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

        // copy back the sigmas
        linker_diagonal_t linker(V.begin());
        // sigmas are the sqrt's of the covariance diagonal
        callLinkFitter<SigmaIdx>(linker, [] (const double& v, double& t) { t = std::sqrt(v); }, vars...);

    }

    /// checkVarSettings

    //std::is_member_function_pointer<decltype(&A::member)>::value

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

    /// getNvars

    template<size_t... Idx, class... Vars_>
    static size_t
    getNvars(const idx_array_t& innerDims, detail::indices<Idx...>, Vars_&&... vars) noexcept {
        using namespace detail;
        return sum_of(
                    std::get<Idx>(innerDims)
                    * // multiply pairwise innerDims by outerDims
                    getOuterDim<Vars_>(
                        std::enable_if<is_stl_cont<Vars>::value>(),
                        std::forward<Vars_>(vars))...
                    );
    }

    /// getOuterDim

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
        return idx_array_t{
            // dispatch depending if Var in parameter pack is container or not
            // if yes, inspect value_type. if no, inspect inner Var itself
            getInnerDim<N, Vars>(std::enable_if<detail::is_stl_cont<Vars>::value>())...
        };
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