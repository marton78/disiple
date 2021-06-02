#pragma once

#include <Eigen/Core>

namespace disiple {

    enum dry_run_t { dry_run };

    template <typename Element,
              typename State, typename Coeffs,
              typename Enable = void>
    class filter_base;

    template <typename Derived, typename State, typename Coeffs>
    class filter_base<Derived, State, Coeffs,
                        typename std::enable_if<
                            !std::is_arithmetic<Derived>::value
                        >::type>
    {
        template <typename> class have_dry_run_apply;

    protected:
        template <typename... Args>
        explicit filter_base(Args&&... args) : coeffs_(std::forward<Args>(args)...) {}

        State& state() { return state_; }
        const State& state() const { return state_; }

        Coeffs& coeffs() { return coeffs_; }
        const Coeffs& coeffs() const { return coeffs_; }

    public:
        /// Initialize filter state to zero.
        void initialize() { state_.initialize(); }

        /// Initialize filter state to a constant (steady state) value.
        /// (As if an infinite number of constant vectors @x_ss@ had
        /// been passed through it.)
        template <typename X>
        void initialize(const Eigen::ArrayBase<X>& x_ss)
        {
            state_.setup(coeffs_, static_cast<int>(x_ss.rows()));
            state_.initialize(coeffs_, x_ss);
        }

        /// Apply filter in-place
        /// @param x Input array, channels-by-time
        template <typename X>
        void apply(Eigen::ArrayBase<X>& x)
        {
            state_.setup(coeffs_, static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                auto xi = x.col(i);
                state_.apply(coeffs_, xi);
            }
        }

        /// Apply filter in-place
        /// @param x Input array, channels-by-time
        template <typename X>
        void apply(Eigen::ArrayBase<X>&& x) { apply(x); }

        /// Apply filter to x and write result to y
        /// @param x Input array, channels-by-time
        /// @param y Output array, channels-by-time
        template <typename X, typename Y>
        void apply(const Eigen::ArrayBase<X>& x, Eigen::ArrayBase<Y>& y)
        {
            state_.setup(coeffs_, static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                auto xi = x.col(i);
                auto yi = y.col(i);
                state_.apply(coeffs_, yi = xi);
            }
        }

        /// Apply filter to x and write result to y
        /// @param x Input array, channels-by-time
        /// @param y Output array, channels-by-time
        template <typename X, typename Y>
        void apply(const Eigen::ArrayBase<X>& x, Eigen::ArrayBase<Y>&& y) { apply(x, y); }

        /// Apply filter, but discard output, just update filter state
        /// @param x Input array, channels-by-time
        template <typename X>
        void apply(const Eigen::ArrayBase<X>& x, dry_run_t)
        {
            state_.setup(coeffs_, static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                enum { tag = have_dry_run_apply<typename X::ColXpr>::value };
                apply_dry(x.col(i), std::integral_constant<bool, tag>());
            }
        }

    private:
        enum { Channels = Eigen::internal::traits<Derived>::RowsAtCompileTime };
        static_assert(Eigen::internal::traits<Derived>::ColsAtCompileTime == 1,
                      "This class is for vectors only");

        template <typename A>
        class have_dry_run_apply
        {
            template <typename S> static uint32_t test(decltype(
                    std::declval<S>().apply(std::declval<Coeffs>(),
                                            std::declval<A const&>(),
                                            dry_run)
                )* = 0);
            template <typename S> static uint16_t test(...);

        public:
            enum { value = sizeof(test<State>(0)) == 4 };
        };

        template <typename X>
        void apply_dry(const Eigen::ArrayBase<X>& xi, std::true_type)
        {
            state_.apply(coeffs_, xi, dry_run);
        }

        template <typename X>
        void apply_dry(const Eigen::ArrayBase<X>& xi, std::false_type)
        {
            state_.apply(coeffs_, backup_ = xi);
        }

        State  state_;
        Coeffs coeffs_;
        Eigen::Array<typename Derived::Scalar, Channels, 1> backup_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };




    template <typename Scalar, typename State, typename Coeffs>
    class filter_base<Scalar, State, Coeffs,
        typename std::enable_if<
            std::is_arithmetic<Scalar>::value
        >::type>
        : public filter_base<Eigen::Array<Scalar, 1, 1>, State, Coeffs>
    {
        typedef Eigen::Array<Scalar, 1, 1>                 array_t;
        typedef Eigen::Map<array_t>                        map_t;
        typedef filter_base<array_t, State, Coeffs>        base_t;

    protected:
        template <typename... Args>
        explicit filter_base(Args&&... args) : base_t(std::forward<Args>(args)...) {}

        using base_t::state;
        using base_t::coeffs;

    public:
        using base_t::initialize;

        Scalar operator()(Scalar x)
        {
            Scalar y;
            base_t::apply(map_t(&x), map_t(&y));
            return y;
        }

        void apply(Scalar x)
        {
            base_t::apply(map_t(&x));
        }

        void apply(Scalar x, Scalar& y)
        {
            base_t::apply(map_t(&x), map_t(&y));
        }

        void apply(Scalar x, dry_run_t)
        {
            base_t::apply(map_t(&x), dry_run);
        }
    };






    template <typename Element, typename Enable = void>
    struct element_traits;

    template <typename T>
    struct element_traits<T, typename std::enable_if<
                !std::is_arithmetic<T>::value
            >::type>
    {
        static_assert(Eigen::internal::traits<T>::ColsAtCompileTime == 1,
                      "This class is for vectors only");
        enum { Channels = Eigen::internal::traits<T>::RowsAtCompileTime };
        typedef typename T::Scalar Scalar;
    };

    template <typename T>
    struct element_traits<T, typename std::enable_if<
                std::is_arithmetic<T>::value
            >::type>
    {
        enum { Channels = 1 };
        typedef T Scalar;
    };


}
