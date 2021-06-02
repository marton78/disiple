#pragma once

#include <disiple/impl/filter_base.hpp>

namespace disiple {

    template <typename, typename, typename> struct band_rms_coeffs;
    template <typename, typename, typename> struct band_rms_state;
    template <typename, typename, typename> struct band_rms_algo;

    /// A filter that calculates the root mean square of a signal
    /// in a given band. The result is corrected via mutliplying
    /// it by sqrt(2) thus giving correct results for sine waves
    /// e.g. band_rms<Array4f, iir<Array4f, 4>, moving_average<Array4f> >
    /// See: http://en.wikipedia.org/wiki/Root_mean_square#RMS_of_common_waveforms

    template <typename Element, typename Bandpass, typename Expectation>
    struct band_rms : public filter_base<Element,
        band_rms_state <typename element_traits<Element>::Scalar, Bandpass, Expectation>,
        band_rms_coeffs<typename element_traits<Element>::Scalar, Bandpass, Expectation>
    >
    {
        enum { Channels = element_traits<Element>::Channels };
        typedef typename element_traits<Element>::Scalar        Scalar;
        typedef band_rms_state <Scalar, Bandpass, Expectation>  state_type;
        typedef band_rms_coeffs<Scalar, Bandpass, Expectation>  coeffs_type;
        typedef filter_base<Element, state_type, coeffs_type>   base_type;

        band_rms() {}

        template <typename TB, typename TE>
        band_rms(TB&& b, TE&& e) : base_type(std::forward<TB>(b), std::forward<TE>(e))
        {}
    };


    template <typename Scalar, typename Bandpass, typename Expectation>
    struct band_rms_coeffs
    {
        band_rms_coeffs() {}

        template <typename TB, typename TE>
        explicit band_rms_coeffs(TB&& b, TE&& e)
        : bpc_(std::forward<TB>(b)), exc_(std::forward<TE>(e)) {}

        Scalar scaling() { return bpc_.scaling() * exc_.scaling(); }

        typename Bandpass::coeffs_type    bpc_;
        typename Expectation::coeffs_type exc_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    template <typename Scalar, typename Bandpass, typename Expectation>
    struct band_rms_state
    {
        typedef band_rms_coeffs<Scalar, Bandpass, Expectation> coeffs_type;

        void setup(const coeffs_type& coeffs, int nchans)
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
        void initialize(const coeffs_type& coeffs, const Eigen::ArrayBase<X>& x_ss)
        {
            bps_.initialize(coeffs.bpc_, x_ss);
            exs_.initialize(coeffs.exc_, x_ss);
        }

        template <typename X>
        void apply(const band_rms_coeffs<Scalar, Bandpass, Expectation>& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            bps_.apply(coeffs.bpc_, xi);
            xi *= xi;
            exs_.apply(coeffs.exc_, xi);
            xi = ( xi.abs() * Scalar(2) ).sqrt();
        }


        typename Bandpass::state_type    bps_;
        typename Expectation::state_type exs_;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

}
