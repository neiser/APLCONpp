#pragma once

#include <type_traits>
#include <tuple>
#include <vector>

namespace APLCON_ {

namespace {

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


// simple index_sequence, as it's not available in C++11
template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};

// compare two std::array contents at compile time
template<class T, std::size_t N, std::size_t... Idx>
constexpr bool compare_array_impl(const std::array<T,N> a, const std::array<T,N> b, indices<Idx...>) noexcept {
    return std::make_tuple(a[Idx]...) == std::make_tuple(b[Idx]...);
}

template<class T, class U, size_t N>
constexpr bool compare_array(const std::array<T,N>& a, const std::array<U,N>& b) noexcept {
    return compare_array_impl(a, b, build_indices<N>());
}

// add just adds the given parameter pack together
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

} // end of anonymous namespace

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

static constexpr size_t c_is_constsize = c_is_single    | c_is_array;
static constexpr size_t c_is_multi     = c_is_container | c_is_array;

static_assert((c_is_constsize & c_is_container) == 0, "Logic error");
static_assert((c_is_multi & c_is_single) == 0, "Logic error");

} // namespace APLCON_