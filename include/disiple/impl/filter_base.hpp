#pragma once

#include <Eigen/Core>

namespace disiple {

    enum DryRun { dry_run };

    template <typename Scalar, typename Derived>
    class FilterBase
    {
        static_assert(std::is_arithmetic<Scalar>::value,
                      "Scalar must be an arithmetic value");
        template <typename, typename> class have_dry_run_apply;
        using Map1 = Eigen::Map<Eigen::Array<Scalar, 1, 1>>;

    public:
        auto&       state()        { return static_cast<Derived*>(this)->state_; }
        const auto& state() const  { return static_cast<Derived const*>(this)->state_; }

        auto&       coeffs()       { return static_cast<Derived*>(this)->coeffs_; }
        const auto& coeffs() const { return static_cast<Derived const*>(this)->coeffs_; }

        /// Initialize filter state to zero.
        void initialize() { state().initialize(); }

        /// Initialize filter state to a constant (steady state) value.
        /// (As if an infinite number of constant vectors @x_ss@ had
        /// been passed through it.)
        template <typename X>
        void initialize(const Eigen::ArrayBase<X>& x_ss)
        {
            state().setup(coeffs(), static_cast<int>(x_ss.rows()));
            state().initialize(coeffs(), x_ss);
        }

        /// Apply filter in-place
        /// @param x Input array, channels-by-time
        template <typename X>
        void apply(Eigen::ArrayBase<X>& x)
        {
            state().setup(coeffs(), static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                auto xi = x.col(i);
                state().apply(coeffs(), xi);
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
            state().setup(coeffs(), static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                auto xi = x.col(i);
                auto yi = y.col(i);
                state().apply(coeffs(), yi = xi);
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
            state().setup(coeffs(), static_cast<int>(x.rows()));
            for (Eigen::DenseIndex i=0; i<x.cols(); ++i)
            {
                enum { tag = have_dry_run_apply<decltype(Derived::coeffs_), typename X::ColXpr>::value };
                apply_dry(x.col(i), std::integral_constant<bool, tag>());
            }
        }

        Scalar operator()(Scalar x)
        {
            Scalar y;
            apply(Map1(&x), Map1(&y));
            return y;
        }

        void apply(Scalar x)
        {
            apply(Map1(&x));
        }

        void apply(Scalar x, Scalar& y)
        {
            apply(Map1(&x), Map1(&y));
        }

        void apply(Scalar x, DryRun)
        {
            apply(Map1(&x), dry_run);
        }

    private:
        template <typename C, typename A>
        class have_dry_run_apply
        {
            template <typename S> static uint32_t test(decltype(
                    std::declval<S>().apply(std::declval<C const&>(),
                                            std::declval<A const&>(),
                                            dry_run)
                )* = 0);
            template <typename S> static uint16_t test(...);

        public:
            enum { value = sizeof(test<decltype(Derived::state_)>(0)) == 4 };
        };

        template <typename X>
        void apply_dry(const Eigen::ArrayBase<X>& xi, std::true_type)
        {
            state().apply(coeffs(), xi, dry_run);
        }

        template <typename X>
        void apply_dry(const Eigen::ArrayBase<X>& xi, std::false_type)
        {
            state().apply(coeffs(), backup_ = xi);
        }

        Eigen::Array<Scalar, /*Derived::Channels*/-1, 1> backup_;
    };

}
