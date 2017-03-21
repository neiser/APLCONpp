#include <iostream>
#include <sstream>
#include <bitset>
#include <iomanip>
#include <vector>
#include <array>

#include <cxxabi.h>

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

template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};

template<class T, size_t N, size_t... Idx>
constexpr bool compare_array_impl(const std::array<T,N> a, const std::array<T,N> b, indices<Idx...>) noexcept {
    return std::tie(a[Idx]...) == std::tie(b[Idx]...);
}

template<class T, class U, size_t N>
constexpr bool compare_array(const std::array<T,N>& a, const std::array<U,N>& b) noexcept {
    return compare_array_impl(a, b, build_indices<N>());
}

//template<typename... T>
//std::tuple<const T&...> ctie(const T&... args)
//{
//    return std::tie( args... );
//}

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

static constexpr size_t c_is_nothing    = 0;
static constexpr size_t c_is_container  = 1;
static constexpr size_t c_is_array      = 2;
static constexpr size_t c_is_arithmetic = 3;

template<typename T>
struct constraint_test {
    static constexpr size_t value = is_stl_container_like<T>::value
                                    ? c_is_container
                                    : (std::is_arithmetic<T>::value
                                       ? c_is_arithmetic
                                       : c_is_nothing);
};

template<typename T, std::size_t N_>
struct constraint_test<std::array<T, N_>> {
    static constexpr size_t value = std::is_arithmetic<T>::value ? c_is_array : c_is_nothing;
    static constexpr size_t N = N_;
};

struct linker_linear_t {

    using it_t = std::vector<double>::iterator;
    it_t it;

    constexpr linker_linear_t(const it_t& it_)
        : it(it_) {}

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), void>::type
    operator()(Func, const std::tuple<Tp...>&)
    { }

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t)
    {
        f(*it++, std::get<I>(t));
        this->operator()<Func, I + 1, Tp...>(f, t);
    }

};

struct linker_diagonal_t {
    using it_t = std::vector<double>::iterator;
    it_t it;
    std::size_t n_row = 0;

    constexpr linker_diagonal_t(it_t it_) : it(it_) {}

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), void>::type
    operator()(Func, const std::tuple<Tp...>&)
    { }

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t)
    {
        it = std::fill_n(it, n_row++, 0); // this skips the row and moves iterator to diagonal element
        f(*it++, std::get<I>(t));
        this->operator()<Func, I + 1, Tp...>(f, t);
    }
};

static constexpr auto ValueIdx = 0;
static constexpr auto SigmaIdx = 1;

template<class ... Types>
struct Fitter {

    using X_t = std::vector<double>;
    X_t X;
    X_t V;
    X_t F;

    static constexpr auto Nt = sizeof...(Types);
    using idx_array_t = std::array<std::size_t, Nt>;

    template<class ... Constraints>
    void DoFit(Types&... types, const Constraints&... constraints) {

        const auto nConstraints = getConstraintDim(types..., constraints...);
        cout << "Constraints: " << nConstraints  << endl;

        constexpr auto innerDim = getInnerDim<ValueIdx>(build_indices<Nt>());
        {
            constexpr auto innerDim_sigmas = getInnerDim<SigmaIdx>(build_indices<Nt>());
            static_assert(compare_array(innerDim, innerDim_sigmas), "The implementation of link_fitter returns unequal number of values and sigmas.");
        }
        const auto nVar = getNvars(innerDim, build_indices<Nt>(), types...);

        cout << "Variables: " << nVar  << endl;

        {
            X.resize(nVar);
            linker_linear_t linker(X.begin());
            // sigmas are the sqrt's of the covariance diagonal
            applyto<ValueIdx>(linker, [] (double& v, const double& t) { v = t; }, types...);
        }

        {
            V.resize((nVar*nVar+nVar)/2);
            linker_diagonal_t linker(V.begin());
            applyto<SigmaIdx>(linker, [] (double& v, const double& t) { v = t*t; }, types...);
        }

//        {
//            F.resize(nConstraints);

//        }

        for(auto x : X)
            cout << x << ' ';
        cout << endl;

        for(auto v : V) {
            cout << v << ' ';
        }
        cout << endl;
    }

    /// applyto

    template<size_t N, class Filler, class Func, class Type, class... Types_>
    static constexpr void
    applyto(Filler&& filler, Func&& f, Type&& type, Types_&&... types) noexcept {
        applyto<N>(std::forward<Filler>(filler), std::forward<Func>(f),
                   std::enable_if<is_stl_container_like<Type>::value>(), std::forward<Type>(type));
        applyto<N>(std::forward<Filler>(filler), std::forward<Func>(f), std::forward<Types_>(types)...);
    }

    template<size_t N, class Filler, class Func>
    static constexpr void
    applyto(Filler&&, Func&&) noexcept {}

    template<size_t N, class Filler, class Func, class Type>
    static constexpr void
    applyto(Filler&& filler, Func&& f, std::enable_if<false>, Type&& type) noexcept {
        std::forward<Filler>(filler)(std::forward<Func>(f), std::forward<Type>(type).template link_fitter<N>());
    }

    template<size_t N, class Filler, class Func, class Type>
    static constexpr void
    applyto(Filler&& filler, Func&& f, std::enable_if<true>, Type&& container) noexcept {
        for(auto& type : std::forward<Type>(container))
            applyto<N>(std::forward<Filler>(filler), std::forward<Func>(f), std::enable_if<false>(), type);
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

    template<class Type>
    static size_t
    getOuterDim(std::enable_if<true>, Type&& type) noexcept {
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
        using f_t = decltype(&Type::template link_fitter<N>);
        using r_t = typename std::result_of<f_t(Type)>::type;
        return std::tuple_size<r_t>::value;

    }

    template<size_t N, class Type>
    static constexpr size_t
    getInnerDim(std::enable_if<true>) noexcept {
        return getInnerDim<N, typename Type::value_type>(std::enable_if<false>());
    }

    /// getConstraintDim

    template<class Constraint, class ... Constraints>
    static constexpr std::size_t
    getConstraintDim(const Types&... types, const Constraint& c, const Constraints& ... c_) noexcept {
        return getConstraintDim(types..., c) + getConstraintDim(types..., c_...);
    }

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<typename std::result_of<Constraint(Types...)>::type>::value == c_is_container, std::size_t>::type
    getConstraintDim(const Types&... types, const Constraint& c) noexcept {
        // no other way than calling it
        return c(types...).size();
    }

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<typename std::result_of<Constraint(Types...)>::type>::value == c_is_array, std::size_t>::type
    getConstraintDim(const Types&..., const Constraint&) noexcept {
        return constraint_test<typename std::result_of<Constraint(Types...)>::type>::N;
    }

    template<class Constraint>
    static constexpr typename std::enable_if<constraint_test<typename std::result_of<Constraint(Types...)>::type>::value == c_is_arithmetic, std::size_t>::type
    getConstraintDim(const Types&..., const Constraint&) noexcept {
        return 1;
    }
};

using ValueSigma_t = std::tuple<double,  double>;

struct A {

    ValueSigma_t E, px, py, pz;

    template<size_t N>
    std::tuple<double&,double&,double&,double&> link_fitter() noexcept {
        cout << "A fitter called N=" << N << endl;
        return std::tie(std::get<N>(E), std::get<N>(px), std::get<N>(py), std::get<N>(pz));
    }
};


struct B {
    ValueSigma_t a;

    template<size_t N>
    std::tuple<double&> link_fitter() noexcept {
        cout << "B fitter called N=" << N << endl;
        //f(std::get<N>(a));
        return std::tie(std::get<N>(a));
    }
};


int main()
{

    A a{{1,9},{2,8},{3,7},{4,6}};

    vector<B> b{
        {{1,1}},
        {{5,2}},
        {{8,3}}
    };

    auto constraint1 = [] (const A& a) {
        cout << "Vector was called: " << std::get<ValueIdx>(a.px) << endl;
        return std::vector<double>{6,7};
    };

    auto constraint2 = [] (const A& a, const vector<B>& b) {
        cout << "Array was called: " << std::get<ValueIdx>(a.E) << std::get<ValueIdx>(b.front().a) << endl;
        return std::array<double, 8>();
    };


    Fitter<A> fitter1;
    cout << ">>>> Fitter1" << endl;
    fitter1.DoFit(a, constraint1);
    cout << endl;

    cout << ">>>> Fitter2" << endl;
    Fitter<A, vector<B>> fitter2;
    fitter2.DoFit(a, b, constraint2);
    cout << endl;
}
