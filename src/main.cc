#include <iostream>
#include <sstream>
#include <bitset>
#include <iomanip>
#include <vector>
#include <array>
#include <tuple>
#include <functional>
#include <cmath>

#include <cxxabi.h>

extern "C" {
#include "wrapper/APLCON.h"
}

using namespace std;

// general helpers

template<typename T>
std::string getTypeAsString() {
    return abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
}

template<typename T>
struct is_stl_container_like
{
    template<typename A>
    static constexpr bool test(
            A * pt,
            const A* cpt = nullptr,
            decltype(pt->begin()) * = nullptr,
            decltype(pt->end()) * = nullptr,
            decltype(cpt->begin()) * = nullptr,
            decltype(cpt->end()) * = nullptr) noexcept
    {
        typedef typename A::iterator iterator;
        typedef typename A::const_iterator const_iterator;
        return  std::is_same<decltype(pt->begin()),iterator>::value &&
                std::is_same<decltype(pt->end()),iterator>::value &&
                std::is_same<decltype(cpt->begin()),const_iterator>::value &&
                std::is_same<decltype(cpt->end()),const_iterator>::value;
    }

    template<typename A>
    static constexpr bool test(...) noexcept {
        return false;
    }

    typedef typename std::decay<T>::type test_type;
    static constexpr bool value = test<test_type>(nullptr);

};

// put this into std namespace, as we target T exactly to
// those STL containers for ostream'ing (this is probably the only occasion
// where it makes sense to pollute std namespace)
// clang 3.9 otherwise complains when compiling test, btw
// (gcc 6.1 is happy though even if it's in namespace ant)
namespace std {

template<class T>
// use SFINAE to restrict this templated operator to STL containers such as vector,list,map,set
typename std::enable_if<is_stl_container_like<T>::value, std::ostream&>::type
operator<< (std::ostream& stream, const T& v)
{
    stream << "[";
    for(auto it = std::begin(v); it != std::end(v); ++it)
    {
        stream << *it;
        if(std::next(it) != std::end(v))
            stream << ";";
    }
    stream << "]";
    return stream;
}

// make std::pair printable for std::map support
template<class U, class V>
std::ostream&
operator<< (std::ostream& stream, const std::pair<U, V>& p) {
    return stream << p.first << "=" << p.second;
}

} // namespace std


template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};

template<class T, size_t N, size_t... Idx>
constexpr bool compare_array_impl(const std::array<T,N> a, const std::array<T,N> b, indices<Idx...>) noexcept {
    return std::make_tuple(a[Idx]...) == std::make_tuple(b[Idx]...);
}

template<class T, class U, size_t N>
constexpr bool compare_array(const std::array<T,N>& a, const std::array<U,N>& b) noexcept {
    return compare_array_impl(a, b, build_indices<N>());
}

template<class T = void>
static constexpr size_t
add() noexcept {
    return 0;
}

template<class T, class... Ts>
static constexpr size_t
add(const T& t, const Ts&... ts) noexcept {
    return t + add(ts...);
}

// start of fitter helpers

static constexpr size_t c_is_nothing     = 0;
static constexpr size_t c_is_container   = (1 << 0);
static constexpr size_t c_is_array       = (1 << 1);
static constexpr size_t c_is_single      = (1 << 2);

template<typename T>
struct constraint_test_impl {
    static constexpr size_t value =
            is_stl_container_like<T>::value ? c_is_container
                                            : (std::is_arithmetic<T>::value
                                               ? c_is_single
                                               : c_is_nothing);

    static constexpr size_t N = std::is_arithmetic<T>::value; // zero for non-constsize T
};

template<typename T, std::size_t N_>
struct constraint_test_impl<std::array<T, N_>> {
    static constexpr size_t value = std::is_arithmetic<T>::value ? c_is_array : c_is_nothing;
    static constexpr size_t N = N_;
};

template<typename T, size_t Mask = c_is_nothing>
struct constraint_test;

template<typename F, typename... Args, size_t Mask>
struct constraint_test<F(Args...), Mask> {
    using r_t = typename std::result_of<F(Args...)>::type;
    static constexpr bool value = static_cast<bool>(constraint_test_impl<r_t>::value & Mask);
    constexpr size_t getN() const noexcept {
        return constraint_test_impl<r_t>::N;
    }
};


struct linker_t {
    using it_t = std::vector<double>::iterator;
    it_t it;

    constexpr linker_t(const it_t& it_)
        : it(it_) {}

    // empty recursion end over index I
    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), void>::type
    operator()(Func, const std::tuple<Tp...>&)
    {}
};

struct linker_linear_t : linker_t {
    using linker_t::linker_t;

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t)
    {
        f(*it++, std::get<I>(t));
        linker_t::operator()<Func, I + 1, Tp...>(f, t);
    }
};

struct linker_diagonal_t : linker_t {
    using linker_t::linker_t;

    std::size_t n_row = 0;

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t)
    {
        // this skips the row and moves iterator to diagonal element
        it = std::fill_n(it, n_row++, 0);
        f(*it++, std::get<I>(t));
        linker_t::operator()<Func, I + 1, Tp...>(f, t);
    }
};



static constexpr auto ValueIdx = 0;
static constexpr auto SigmaIdx = 1;

static constexpr size_t c_is_constsize = c_is_single    | c_is_array;
static constexpr size_t c_is_multi     = c_is_container | c_is_array;

static_assert((c_is_constsize & c_is_container) == 0, "Logic error");
static_assert((c_is_multi & c_is_single) == 0, "Logic error");

template<class ... Types>
struct Fitter {

    std::vector<double> X;
    std::vector<double> V;
    std::vector<double> F_dynamic; // only used if non-constexpr constraints size


    static constexpr auto Nt = sizeof...(Types);
    using idx_array_t = std::array<std::size_t, Nt>;

    template<class ... Constraints>
    void DoFit(Types&... types, Constraints&&... constraints) {

        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nt>());
        {
            constexpr auto innerDim_sigmas = getInnerDim<SigmaIdx>(build_indices<Nt>());
            static_assert(compare_array(innerDim, innerDim_sigmas), "The implementation of link_fitter returns unequal number of values and sigmas.");
        }
        const auto nVar = getNvars(innerDim, build_indices<Nt>(), types...);

        {
            X.resize(nVar);
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            callLinkFitter<ValueIdx>(linker, [] (double& v, const double& t) { v = t; }, types...);
        }

        {
            V.resize((nVar*nVar+nVar)/2);
            linker_diagonal_t linker(V.begin());
            callLinkFitter<SigmaIdx>(linker, [] (double& v, const double& t) { v = t*t; }, types...);
        }

        constexpr auto isConstexpr = add(isConstraintsSizeConstexpr<Constraints>()...)==0;

        // dispatch via AllocF
        AllocF(std::enable_if<isConstexpr>(), types..., std::forward<Constraints>(constraints)...);

    }

    template<class ... Constraints>
    void AllocF(std::enable_if<true>, Types&... types, Constraints&&... constraints) {
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
    static typename std::enable_if<constraint_test<Constraint(Types...), c_is_multi>::value, void>::type
    callConstraint(It& it, Constraint&& constraint, const Types&... types) noexcept {
        // Making F_ constexpr is probably too limiting for constraints (think of TreeFitter with rather complex lambdas)
        const auto& F_ = std::forward<Constraint>(constraint)(types...);
        it = std::copy(F_.begin(), F_.end(), it);
    }

    template<class It, class Constraint>
    static typename std::enable_if<constraint_test<Constraint(Types...), c_is_single>::value, void>::type
    callConstraint(It& it, Constraint&& constraint, const Types&... types) noexcept {
        // single scalar is just copied to current position
        *it++ = std::forward<Constraint>(constraint)(types...);
    }

    /// callLinkFitter

    template<size_t N, class Linker, class Func, class Type, class... Types_>
    static void callLinkFitter(Linker&& linker, Func&& f, Type&& type, Types_&&... types) noexcept {
        callLinkFitter<N>(std::forward<Linker>(linker), std::forward<Func>(f),
                          std::enable_if<is_stl_container_like<Type>::value>(), std::forward<Type>(type));
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
    getNvars(const idx_array_t& innerDims, indices<Idx...>, Types_&&... types) noexcept {
        return add(
                    get<Idx>(innerDims)
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
    getInnerDim(indices<Idx...>) noexcept {
        return idx_array_t{getInnerDim<N, Types>(std::enable_if<is_stl_container_like<Types>::value>())...};
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
    static constexpr typename std::enable_if<constraint_test<Constraint(Types...), c_is_constsize>::value, std::size_t>::type
    getConstraintDim(Constraint&&, const Types&...) noexcept {
        return constraint_test<Constraint(Types...)>().getN();
    }

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<Constraint(Types...), c_is_container>::value, std::size_t>::type
    getConstraintDim(Constraint&& constraint, const Types&... types) noexcept {
        // call the constraint with the given types
        return std::forward<Constraint>(constraint)(types...).size();
    }

    /// isConstraintsSizeConstexpr

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<Constraint(Types...), c_is_constsize>::value, std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 0;
    }

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<Constraint(Types...), c_is_container>::value, std::size_t>::type
    isConstraintsSizeConstexpr() noexcept {
        return 1;
    }
};

struct X {

    constexpr X(double v, double s) : Value(v), Sigma(s) {}

    double Value;
    double Sigma;

    template<size_t N>
    std::tuple<double&> linkFitter() noexcept {
        // can this be written easier?
        return std::tie(std::get<N>(std::tie(Value, Sigma)));
    }

    friend std::ostream& operator<<(std::ostream& s, const X& o) {
        return s << "(" << o.Value << "," << o.Sigma << ")";
    }
};

#include "benchmark/benchmark_api.h"

static void BM_ErrorPropagation(benchmark::State& state) {

    auto a_and_b_is_c = [] (const X& a, const X& b, const X& c) {
        return c.Value - a.Value - b.Value;
    };

    while (state.KeepRunning()) {

        X a{10, 0.3};
        X b{20, 0.4};
        X c{ 0, 0.0};

        Fitter<X, X, X> fitter;
        fitter.DoFit(a, b, c, a_and_b_is_c);

    }
}
BENCHMARK(BM_ErrorPropagation);

BENCHMARK_MAIN()

