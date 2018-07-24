#pragma once

#include <type_traits>
#include <tuple>
#include <vector>
#include <cmath>

namespace APLCON {
namespace detail {

namespace {

template<typename T>
struct is_stl_cont
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

template<class T, typename Enable = void>
struct decay_stl_cont {
    using type = T;
};

template<class T>
struct decay_stl_cont<T, typename std::enable_if<is_stl_cont<T>::value>::type> {
    using type = typename T::value_type;
};


template<typename T, typename Ret>
struct has_getFitterSettings
{
    template<typename U, Ret (U::*)(size_t) const> struct SFINAE {};
    template<typename U> static char Test(SFINAE<U, &U::template getFitterSettings<0>>*);
    template<typename U> static int Test(...);
    static constexpr bool value = sizeof(Test<T>(0)) == sizeof(char);
};

template < typename... T >
constexpr std::tuple<const T&...> ctie( const T&... args )
{
    return std::tie( args... );
}

// simple index_sequence, as it's not available in C++11
template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};

// compare two std::array contents at compile time
template<class T, std::size_t N, std::size_t... Idx>
constexpr bool
compare_array_impl(const std::array<T,N> a, const std::array<T,N> b, indices<Idx...>) noexcept {
    return std::make_tuple(a[Idx]...) == std::make_tuple(b[Idx]...);
}

template<class T, class U, std::size_t N>
constexpr bool
compare_array(const std::array<T,N>& a, const std::array<U,N>& b) noexcept {
    return compare_array_impl(a, b, build_indices<N>());
}

// sum_of just adds the given parameter pack together
template<class T = void>
static constexpr std::size_t
sum_of() noexcept {
    return 0;
}

template<class T, class... Ts>
static constexpr std::size_t
sum_of(const T& t, const Ts&... ts) noexcept {
    return static_cast<size_t>(t) + sum_of(ts...);
}

template<class T, std::size_t N, std::size_t... Idx>
static constexpr std::size_t
sum_of_array_impl(const std::array<T,N>& arr, indices<Idx...>) noexcept {
    return sum_of(std::get<Idx>(arr)...);
}

template<class T, std::size_t N>
static constexpr std::size_t
sum_of_array(const std::array<T,N>& arr) noexcept {
    return sum_of_array_impl(arr, build_indices<N>());
}

template<class T, std::size_t N, std::size_t... Idx>
static constexpr std::array<T, sizeof...(Idx)>
prod_of_array_impl(const std::array<T,N>& a, const std::array<T,N>& b, indices<Idx...>) noexcept {
    return {{std::get<Idx>(a)*std::get<Idx>(b)...}};
}

template<class T, std::size_t N>
static constexpr std::array<T,N>
prod_of_array(const std::array<T,N>& a, const std::array<T,N>& b) noexcept {
    return prod_of_array_impl(a, b, build_indices<N>());
}

template<typename Arg = void>
bool check_args() {
    return true;
}

template<typename Arg, typename... Args>
bool check_args(Arg&& arg, Args&&... args) {
    if(std::is_convertible<Arg, double>::value)
        return std::isfinite(std::forward<Arg>(arg)) && check_args(std::forward<Args>(args)...);
    if(std::is_convertible<Arg, int>::value)
        return std::forward<Arg>(arg) >= 0 && check_args(std::forward<Args>(args)...);
    return false;
}

template<typename Func, class... Args>
void call_if_set(Func&& func, Args&&... args)  noexcept {
    if(check_args(std::forward<Args>(args)...))
        std::forward<Func>(func)(args...);
}


} // end of anonymous namespace

static constexpr std::size_t c_is_nothing     = 0;
static constexpr std::size_t c_is_container   = (1 << 0);
static constexpr std::size_t c_is_array       = (1 << 1);
static constexpr std::size_t c_is_single      = (1 << 2);

template<typename T>
struct constraint_test_impl {
    static constexpr std::size_t value =
            is_stl_cont<T>::value ? c_is_container
                                            : (std::is_arithmetic<T>::value
                                               ? c_is_single
                                               : c_is_nothing);

    static constexpr std::size_t N = std::is_arithmetic<T>::value; // zero for non-constsize T
    // give some better error message in case of non-usable constraint
    static_assert(value != c_is_nothing, "You probably passed a constraint with did not return a STL container or arithmetic type");
};

template<typename T, std::size_t N_>
struct constraint_test_impl<std::array<T, N_>> {
    static constexpr std::size_t value = std::is_arithmetic<T>::value ? c_is_array : c_is_nothing;
    static constexpr std::size_t N = N_;
};

template<typename T, std::size_t Mask = c_is_nothing>
struct constraint_test {
    static_assert(std::integral_constant<T, false>::value,
                  "First template parameter T needs to be of function type.");
};

template<typename F, typename... Args, std::size_t Mask>
struct constraint_test<F(Args...), Mask> {
    using r_t = typename std::result_of<F(Args...)>::type;
    static constexpr bool value = static_cast<bool>(constraint_test_impl<r_t>::value & Mask);
    constexpr std::size_t getN() const noexcept {
        return constraint_test_impl<r_t>::N;
    }
};

struct linker_t {
    using it_t = std::vector<double>::iterator;
    it_t it;
    constexpr linker_t(const it_t& it_) : it(it_) {}
};

struct linker_linear_t : linker_t {
    using linker_t::linker_t;

    // empty recursion end over index I
    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), void>::type
    operator()(Func, const std::tuple<Tp...>&) noexcept
    {}

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t) noexcept
    {
        f(*it++, std::get<I>(t));
        operator()<Func, I + 1, Tp...>(f, t);
    }
};

struct linker_diagonal_t : linker_t {
    using linker_t::linker_t;

    std::size_t n_row = 0;

    // empty recursion end over index I
    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I == sizeof...(Tp), void>::type
    operator()(Func, const std::tuple<Tp...>&) noexcept
    {}

    template<typename Func, std::size_t I = 0, typename... Tp>
    typename std::enable_if<I < sizeof...(Tp), void>::type
    operator()(Func f, const std::tuple<Tp...>& t) noexcept
    {
        // this skips the row and moves iterator to diagonal element
        it = std::fill_n(it, n_row++, 0);
        f(*it++, std::get<I>(t));
        operator()<Func, I + 1, Tp...>(f, t);
    }
};

static constexpr std::size_t c_is_constsize = c_is_single    | c_is_array;
static constexpr std::size_t c_is_multi     = c_is_container | c_is_array;

static_assert((c_is_constsize & c_is_container) == 0, "Logic error");
static_assert((c_is_multi & c_is_single) == 0, "Logic error");

}} // namespace APLCON::detail