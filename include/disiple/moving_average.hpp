#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/mavg_impl.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <typename Scalar, typename... Options>
    class MovingAverage : public FilterBase<Scalar, MovingAverage<Scalar, Options...>>,
                          public Parameters<
                                List<Options...>,
                                OptionalValue<int, Length, Eigen::Dynamic>,
                                OptionalValue<int, Stages, 1>,
                                OptionalValue<int, Channels, 1>
                            >
    {
    public:
        using State  = MAvgState <Scalar, MovingAverage::LengthValue::Static, MovingAverage::ChannelsValue::Static, MovingAverage::StagesValue::Static>;
        using Coeffs = MAvgCoeffs<Scalar, MovingAverage::LengthValue::Static, MovingAverage::StagesValue::Static>;

        MovingAverage() {}
        explicit MovingAverage(int length) : coeffs_(length) {}
        MovingAverage(int length, int stages) : coeffs_(length, stages) {}

    private:
        friend FilterBase<Scalar, MovingAverage<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };


    template <typename Scalar, typename... Options>
    class CumMovingAverage : public FilterBase<Scalar, CumMovingAverage<Scalar, Options...>>,
                             public Parameters<
                                    List<Options...>,
                                    OptionalValue<int, Channels, 1>
                                >
    {
    public:
        using State  = CumMAvgState<Scalar, CumMovingAverage::ChannelsValue::Static>;
        using Coeffs = CumMAvgCoeffs<Scalar>;

        CumMovingAverage() {}
        explicit CumMovingAverage(int length) : coeffs_(length) {}

    private:
        friend FilterBase<Scalar, CumMovingAverage<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

}
