#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <typename, typename, typename> struct BandRMSCoeffs;
    template <typename, typename, typename> struct BandRMSState;

    /// A filter that calculates the root mean square of a signal
    /// in a given band. The result is corrected via mutliplying
    /// it by sqrt(2) thus giving correct results for sine waves
    /// e.g. BandRMS<Array4f, iir<Array4f, 4>, moving_average<Array4f> >
    /// See: http://en.wikipedia.org/wiki/Root_mean_square#RMS_of_common_waveforms

    template <typename... Options>
    using BandRMSParameters = Parameters<
        List<Options...>,
        OptionalValue<int, Channels, 1>
    >;

    template <typename Scalar, typename Bandpass, typename Expectation, typename... Options>
    struct BandRMS : public FilterBase<Scalar,
        BandRMSParameters<Options...>::channels,
        BandRMSState <Scalar, Bandpass, Expectation>,
        BandRMSCoeffs<Scalar, Bandpass, Expectation>
    >
    {
        using P      = BandRMSParameters<Options...>;
        using State  = BandRMSState <Scalar, Bandpass, Expectation>;
        using Coeffs = BandRMSCoeffs<Scalar, Bandpass, Expectation>;
        using Base   = FilterBase<Scalar, P::channels, State, Coeffs>;

        BandRMS() {}

        template <typename TB, typename TE>
        BandRMS(TB&& b, TE&& e) : Base(std::forward<TB>(b), std::forward<TE>(e))
        {}
    };


    template <typename Scalar, typename Bandpass, typename Expectation>
    struct BandRMSCoeffs
    {
        BandRMSCoeffs() {}

        template <typename TB, typename TE>
        explicit BandRMSCoeffs(TB&& b, TE&& e)
        : bpc_(std::forward<TB>(b)), exc_(std::forward<TE>(e)) {}

        Scalar scaling() { return bpc_.scaling() * exc_.scaling(); }

        typename Bandpass::Coeffs    bpc_;
        typename Expectation::Coeffs exc_;
    };

    template <typename Scalar, typename Bandpass, typename Expectation>
    struct BandRMSState
    {
        using Coeffs = BandRMSCoeffs<Scalar, Bandpass, Expectation>;

        void setup(const Coeffs& coeffs, int nchans)
        {
            bps_.setup(coeffs.bpc_, nchans);
            exs_.setup(coeffs.exc_, nchans);
        }

        void initialize()
        {
            bps_.initialize();
            exs_.initialize();
        }

        template <typename X>
        void initialize(const Coeffs& coeffs, const Eigen::ArrayBase<X>& x_ss)
        {
            bps_.initialize(coeffs.bpc_, x_ss);
            exs_.initialize(coeffs.exc_, x_ss);
        }

        template <typename X>
        void apply(const BandRMSCoeffs<Scalar, Bandpass, Expectation>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            bps_.apply(coeffs.bpc_, xi);
            xi *= xi;
            exs_.apply(coeffs.exc_, xi);
            xi = ( xi.abs() * Scalar(2) ).sqrt();
        }


        typename Bandpass::State    bps_;
        typename Expectation::State exs_;
    };

}
