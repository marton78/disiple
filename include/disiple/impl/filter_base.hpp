#pragma once

#include <Eigen/Core>

namespace disiple {

    enum DryRun { dry_run };

    template <typename Scalar, int Channels,
              typename State, typename Coeffs>
    class FilterBase
    {
        static_assert(std::is_arithmetic<Scalar>::value,
                      "Scalar must be an arithmetic value");
        template <typename> class have_dry_run_apply;
        using Map1 = Eigen::Map<Eigen::Array<Scalar, 1, 1>>;

    protected:
        template <typename... Args>
        explicit FilterBase(Args&&... args) : coeffs_(std::forward<Args>(args)...) {}

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
        void apply(const Eigen::ArrayBase<X>& x, DryRun)
        {
            state_.setup(coeffs_, static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                enum { tag = have_dry_run_apply<typename X::ColXpr>::value };
                apply_dry(x.col(i), std::integral_constant<bool, tag>());
            }
        }

        template <int C = Channels, typename Enable = typename std::enable_if<C == 1>::type>
        Scalar operator()(Scalar x, Enable* = 0)
        {
            Scalar y;
            apply(Map1(&x), Map1(&y));
            return y;
        }

        template <int C = Channels, typename Enable = typename std::enable_if<C == 1>::type>
        void apply(Scalar x, Enable* = 0)
        {
            apply(Map1(&x));
        }

        template <int C = Channels, typename Enable = typename std::enable_if<C == 1>::type>
        void apply(Scalar x, Scalar& y, Enable* = 0)
        {
            apply(Map1(&x), Map1(&y));
        }

        template <int C = Channels, typename Enable = typename std::enable_if<C == 1>::type>
        void apply(Scalar x, DryRun, Enable* = 0)
        {
            apply(Map1(&x), dry_run);
        }

    private:
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
        Eigen::Array<Scalar, Channels, 1> backup_;
    };

}
