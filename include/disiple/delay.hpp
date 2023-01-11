#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/delay_impl.hpp>
#include <disiple/named_params.hpp>


namespace disiple {

    template <typename Scalar, typename... Options>
    class Delay : public FilterBase<Scalar, Delay<Scalar, Options...>>,
                  public Parameters<
                        List<Options...>,
                        OptionalValue<int, Length, Eigen::Dynamic>,
                        OptionalValue<int, Channels, 1>
                    >
    {
    public:
        using State  = DelayState<Scalar, Delay::length, Delay::channels>;
        using Coeffs = DelayCoeffs<Scalar, Delay::length>;

        Delay() {}
        explicit Delay(int length) : coeffs_(length) {}

    private:
        friend FilterBase<Scalar, Delay<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

}
