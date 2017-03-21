#include <iostream>
#include <sstream>
#include <bitset>
#include <iomanip>
#include <vector>
#include <array>

#include <cxxabi.h>

using namespace std;

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
            decltype(cpt->end()) * = nullptr)
    {
        typedef typename A::iterator iterator;
        typedef typename A::const_iterator const_iterator;
        return  std::is_same<decltype(pt->begin()),iterator>::value &&
                std::is_same<decltype(pt->end()),iterator>::value &&
                std::is_same<decltype(cpt->begin()),const_iterator>::value &&
                std::is_same<decltype(cpt->end()),const_iterator>::value;
    }

    template<typename A>
    static constexpr bool test(...) {
        return false;
    }

    typedef typename std::remove_const<T>::type test_type;
    static constexpr bool value = test<test_type>(nullptr);

};


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

struct get_num_args_t {
    std::size_t num_args = 0;
    template<class... Args>
    void operator()(Args&&...) {
        num_args = sizeof...(Args);
    }
};

struct fill_linear_t {

    using it_t = std::vector<double>::iterator;
    it_t x;

    constexpr fill_linear_t(const it_t& x_) : x(x_) {}

    template<class End = void>
    constexpr void operator()() const noexcept {}

    template<class Head, class... Tail>
    void operator()(const Head& head, const Tail&... tail) noexcept {
        *x++ = head;
        this->operator()(tail...);
    }
};

struct fill_diagonal_t {
    using it_t = std::vector<double>::iterator;
    it_t v;
    std::size_t n = 0;

    constexpr fill_diagonal_t(it_t v_) : v(v_) {}

    template<class End = void>
    constexpr void operator()() const noexcept {}

    template<class Head, class... Tail>
    void operator()(const Head& head, const Tail&... tail) noexcept {
        for(std::size_t i = 0; i<n ; i++)
            *v++ = 0;
        *v++ = head;
        n++;
        this->operator()(tail...);
    }
};

template<class ... Types>
struct Fitter {

    using X_t = std::vector<double>;
    X_t X;
    X_t V;
    X_t F;

    template<class ... Constraints>
    void DoFit(Types&... types, const Constraints&... constraints) {
        const auto nConstraints = getConstraintDim(types..., constraints...);
        cout << "Constraints: " << nConstraints  << endl;

        const auto nVar = getTypesDim(types...);
        cout << "Variables: " << nVar  << endl;

        {
            X.resize(nVar);
            fill_linear_t f(X.begin());
            fillX(f, types...);
        }

        {
            V.resize((nVar*nVar+nVar)/2);
            fill_diagonal_t f(V.begin());
            fillV(f, types...);
        }

        {
            F.resize(nConstraints);

        }

        for(auto x : X)
            cout << x << ' ';
        cout << endl;

        for(auto v : V) {
            cout << v << ' ';
        }
        cout << endl;
    }

    /// fillF




    /// fillV

    template<class Filler, class Type, class... Types_>
    static constexpr void
    fillV(Filler& f, const Type& type, const Types_&... types) noexcept {
        fillV(f, std::enable_if<is_stl_container_like<Type>::value>(), type);
        fillV(f, types...);
    }

    template<class Filler>
    static constexpr void
    fillV(Filler&) noexcept {}

    template<class Filler, class Type>
    static constexpr void
    fillV(Filler& f, std::enable_if<false>, const Type& type) noexcept {
        const_cast<Type&>(type).fitter_sigmas(f); // difference from X to V
    }

    template<class Filler, class Type>
    static constexpr typename std::enable_if<is_stl_container_like<Type>::value, void>::type
    fillV(Filler& f, std::enable_if<true>, const Type& container) noexcept {
        for(const auto& type : container)
            fillV(f, std::enable_if<false>(), type);
    }

    /// fillX

    template<class Filler, class Type, class... Types_>
    static constexpr void
    fillX(Filler& f, const Type& type, const Types_&... types) noexcept {
        fillX(f, std::enable_if<is_stl_container_like<Type>::value>(), type);
        fillX(f, types...);
    }

    template<class Filler>
    static constexpr void
    fillX(Filler&) noexcept {}

    template<class Filler, class Type>
    static constexpr void
    fillX(Filler& f, std::enable_if<false>, const Type& type) noexcept {
        const_cast<Type&>(type).fitter_values(f); // difference from X to V
    }

    template<class Filler, class Type>
    static constexpr typename std::enable_if<is_stl_container_like<Type>::value, void>::type
    fillX(Filler& f, std::enable_if<true>, const Type& container) noexcept {
        for(const auto& type : container)
            fillX(f, std::enable_if<false>(), type);
    }

    /// getTypesDim

    template<class Type, class... Types_>
    static constexpr std::size_t
    getTypesDim(const Type& type, const Types_&... types) noexcept {
        return getTypesDim(std::enable_if<is_stl_container_like<Type>::value>(), type) + getTypesDim(types...);
    }

    static constexpr std::size_t
    getTypesDim() noexcept {
        return 0;
    }

    template<class Type>
    static constexpr std::size_t
    getTypesDim(std::enable_if<false>, const Type& type) noexcept {
        // use lambda to have function C++11 constexpr compliant
        return [] (Type& t) {
            get_num_args_t a;
            t.fitter_values(a);
            return a.num_args;
        }(const_cast<Type&>(type));
    }

    template<class Type>
    static constexpr typename std::enable_if<is_stl_container_like<Type>::value, size_t>::type
    getTypesDim(std::enable_if<true>, const Type& container) noexcept {
        // don't know if checking for empty here is really a good idea...
        return container.empty() ? 0 : container.size()*getTypesDim(std::enable_if<false>(), container.front());
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

struct ValueSigma_t {
    double V;
    double S;
};

struct A {

    ValueSigma_t E, px, py, pz;

    template<class Fitter>
    void fitter_values(Fitter& f) {
        cout << "A fitter_values called" << endl;
        f(E.V, px.V, py.V, pz.V);
    }

    template<class Fitter>
    constexpr void fitter_sigmas(Fitter& f) {
        cout << "A fitter_sigmas called" << endl;
        f(E.S, px.S, py.S, pz.S);
    }
};


struct B {
    ValueSigma_t a;

    template<class Fitter>
    void fitter_values(Fitter& f) {
        cout << "B fitter_values called" << endl;
        f(a.V);
    }

    template<class Fitter>
    void fitter_sigmas(Fitter& f) {
        cout << "B fitter_sigmas called" << endl;
        f(a.S);
    }
};


int main()
{

    A a{{0,9},{0,8},{0,7},{0,6}};

    vector<B> b{
        {{1,1}},
        {{5,2}},
        {{8,3}}
    };

    auto constraint1 = [] (const A& a) {
        cout << "Vector was called: " << a.px.S << endl;
        return std::vector<double>{6,7};
    };

    auto constraint2 = [] (const A& a, const vector<B>& b) {
        cout << "Array was called: " << a.E.V << b.front().a.V << endl;
        return std::array<double, 10>();
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
