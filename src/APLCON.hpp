#pragma once

// link to the FORTRAN world
extern "C" {
#include "wrapper/APLCON.h"
}

// some template in extra header file
#include "detail/APLCON_hpp.hpp"
#include <cmath> // for sqrt

template<class ... Types>
struct APLCON {

    // use those constants in your classes to
    // identify what is request in the call to "linkFitter"
    static constexpr auto ValueIdx = 0;
    static constexpr auto SigmaIdx = 1;

    template<class ... Constraints>
    void DoFit(Types&... types, Constraints&&... constraints) {

        using namespace APLCON_;

        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nt>());
        {
            constexpr auto innerDim_sigmas = getInnerDim<SigmaIdx>(build_indices<Nt>());
            static_assert(compare_array(innerDim, innerDim_sigmas), "The implementation of link_fitter returns unequal number of values and sigmas.");
        }
        const auto nVar = getNvars(innerDim, build_indices<Nt>(), types...);

        // fill the initial values
        {
            X.resize(nVar);
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (double& v, const double& t) { v = t; }, types...);
        }

        // fill the initial uncertainties / sigmas
        {
            V.resize((nVar*nVar+nVar)/2);
            linker_diagonal_t linker(V.begin());
            callLinkFitter<SigmaIdx>(linker, [] (double& v, const double& t) { v = t*t; }, types...);
        }


        // dispatch via AllocF to ensure optimal constraint storage handling
        // (and then reach RunAPLCON)
        constexpr auto isConstexpr = add(isConstraintsSizeConstexpr<Constraints>()...)==0;
        AllocF(std::enable_if<isConstexpr>(), types..., std::forward<Constraints>(constraints)...);

    }

private:

    // keep storage here to avoid re-allocations
    std::vector<double> X;
    std::vector<double> V;
    std::vector<double> F_dynamic; // only used if non-constexpr constraints size

    static constexpr auto Nt = sizeof...(Types);
    using idx_array_t = std::array<std::size_t, Nt>;

    template<class ... Constraints>
    void AllocF(std::enable_if<true>, Types&... types, Constraints&&... constraints) {
        using namespace APLCON_;
        constexpr auto nConstraints = add(constraint_test<Constraints(Types...)>().getN()...);
        // alloc on stack
        std::array<double, nConstraints> F_static;
        RunAPLCON(F_static, types..., std::forward<Constraints>(constraints)...);
    }

    template<class ... Constraints>
    void AllocF(std::enable_if<false>, Types&... types, Constraints&&... constraints) {
        const auto nConstraints = add(getConstraintDim<Constraints>(std::forward<Constraints>(constraints), types...)...);
        // alloc on heap
        F_dynamic.resize(nConstraints);
        RunAPLCON(F_dynamic, types..., std::forward<Constraints>(constraints)...);
    }

    template<class F_t, class ... Constraints>
    void RunAPLCON(F_t& F, Types&... types, Constraints&&... constraints) {
        using namespace APLCON_;

        c_aplcon_aplcon(X.size(), F.size());

        c_aplcon_aprint(6, 0); // default output on LUNP 6 (STDOUT)

        // the main convergence loop
        int aplcon_ret = -1;
        do {

            // compute F from constraints
            callConstraints(F.begin(), types..., std::forward<Constraints>(constraints)...);

            // call APLCON iteration
            c_aplcon_aploop(X.data(), V.data(), F.data(), &aplcon_ret);

            // copy X into types, for constraint evaluation
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (const double& v, double& t) { t = v; }, types...);

        }
        while(aplcon_ret<0);

        // copy back the sigmas
        linker_diagonal_t linker(V.begin());
        // sigmas are the sqrt's of the covariance diagonal
        callLinkFitter<SigmaIdx>(linker, [] (const double& v, double& t) { t = std::sqrt(v); }, types...);

    }

    /// callConstraints

    template<class It, class Constraint, class... Constraints>
    static void callConstraints(It it, const Types&... types, Constraint&& constraint, Constraints&&... constraints) noexcept {
        callConstraint<It>(it, std::forward<Constraint>(constraint), types...);
        callConstraints(it, types..., std::forward<Constraints>(constraints)...);
    }

    template<class It>
    static void callConstraints(It, const Types&...)  noexcept {}

    template<class It, class Constraint>
    static typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_multi>::value,
    void>::type
    callConstraint(It& it, Constraint&& constraint, const Types&... types) noexcept {
        // Making F_ constexpr is probably too limiting for constraints (think of TreeFitter with rather complex lambdas)
        const auto& F_ = std::forward<Constraint>(constraint)(types...);
        it = std::copy(F_.begin(), F_.end(), it);
    }

    template<class It, class Constraint>
    static typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_single>::value,
    void>::type
    callConstraint(It& it, Constraint&& constraint, const Types&... types) noexcept {
        // single scalar is just copied to current position
        *it++ = std::forward<Constraint>(constraint)(types...);
    }

    /// callLinkFitter

    template<size_t N, class Linker, class Func, class Type, class... Types_>
    static void callLinkFitter(Linker&& linker, Func&& f, Type&& type, Types_&&... types) noexcept {
        callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f),
                          std::enable_if<APLCON_::is_stl_container_like<Type>::value>(), std::forward<Type>(type));
        callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f), std::forward<Types_>(types)...);
    }

    template<size_t N, class Linker, class Func>
    static void callLinkFitter(Linker&&, Func&&) noexcept {}

    template<size_t N, class Linker, class Func, class Type>
    static void callLinkFitter(Linker&& linker, Func&& f, std::enable_if<false>, Type&& type) noexcept {
        std::forward<Linker>(linker)(std::forward<Func>(f), std::forward<Type>(type).template linkFitter<N>());
    }

    template<size_t N, class Linker, class Func, class Type>
    static void callLinkFitter(Linker&& linker, Func&& f, std::enable_if<true>, Type&& container) noexcept {
        for(auto& type : std::forward<Type>(container))
            callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f), std::enable_if<false>(), type);
    }

    /// getNvars

    template<size_t... Idx, class... Types_>
    static size_t
    getNvars(const idx_array_t& innerDims, APLCON_::indices<Idx...>, Types_&&... types) noexcept {
        using namespace APLCON_;
        return add(
                    std::get<Idx>(innerDims)
                    * // multiply pairwise innerDims by outerDims
                    getOuterDim<Types_>(
                        std::enable_if<is_stl_container_like<Types>::value>(),
                        std::forward<Types_>(types))...
                    );
    }

    /// getOuterDim

    template<class Type>
    static size_t
    getOuterDim(std::enable_if<true>, Type&& type) noexcept {
        // no other way than calling non-constexpr .size()
        return std::forward<Type>(type).size();
    }

    template<class Type>
    static constexpr size_t
    getOuterDim(std::enable_if<false>, Type&&) noexcept {
        return 1;
    }

    /// getInnerDim

    template<size_t N, size_t... Idx>
    static constexpr idx_array_t
    getInnerDim(APLCON_::indices<Idx...>) noexcept {
        return idx_array_t{getInnerDim<N, Types>(std::enable_if<APLCON_::is_stl_container_like<Types>::value>())...};
    }

    template<size_t N, class Type>
    static constexpr std::size_t
    getInnerDim(std::enable_if<false>) noexcept {
        // inspect the templated function
        using f_t = decltype(&Type::template linkFitter<N>);
        using r_t = typename std::result_of<f_t(Type)>::type;
        return std::tuple_size<r_t>::value;

    }

    template<size_t N, class Type>
    static constexpr size_t
    getInnerDim(std::enable_if<true>) noexcept {
        return getInnerDim<N, typename Type::value_type>(std::enable_if<false>());
    }


    /// getConstraintDim

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_constsize>::value,
    std::size_t>::type
    getConstraintDim(Constraint&&, const Types&...) noexcept {
        return APLCON_::constraint_test<Constraint(Types...)>().getN();
    }

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_container>::value,
    std::size_t>::type
    getConstraintDim(Constraint&& constraint, const Types&... types) noexcept {
        // call the constraint with the given types
        return std::forward<Constraint>(constraint)(types...).size();
    }

    /// isConstraintsSizeConstexpr

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_constsize>::value,
    std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 0;
    }

    template<class Constraint>
    static constexpr typename
    std::enable_if<
    APLCON_::constraint_test<Constraint(Types...), APLCON_::c_is_container>::value,
    std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 1;
    }
};