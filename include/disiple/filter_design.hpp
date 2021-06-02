#pragma once

#include <disiple/impl/pole_zero_pair.hpp>
#include <utility>
#include <vector>
#include <Eigen/Core>

namespace disiple {

    //
    // Describes an analog filter as a collection of poles and zeros along with
    // normalization information to achieve a specified gain at a specified
    // frequency. The poles and zeros lie in the s plane.
    //

    struct iir_prototype;

    //
    // Describes a digital filter as a collection of poles and zeros along with
    // normalization information to achieve a specified gain at a specified
    // frequency. The poles and zeros lie in the z plane.
    //

    class iir_design
    {
    public:
        iir_design() : num_poles_(0) { pairs_.reserve(16); }

        template <typename Transform>
        iir_design(const iir_prototype& p, Transform t);

        void add_single(double pole, double zero);
        void add_conjugate_pair(const complex_t pole, const complex_t zero);
        void add(const complex_pair& poles, const complex_pair& zeros);
        void add_notch(double Wo, double BW);

        int num_poles()                          const { return num_poles_; }
        const pole_zero_pair& operator[] (int i) const { return pairs_[i]; }
        double normal_w()                        const { return normal_w_; }
        double normal_gain()                     const { return normal_g_; }
        void set_normal(double w, double g) { normal_w_ = w; normal_g_ = g; }

    private:
        int num_poles_;
        std::vector<pole_zero_pair> pairs_;
        double normal_w_;
        double normal_g_;
    };

    iir_design notch(double Wo, double BW);
    iir_design comb(int N, double BW);

    struct iir_prototype : iir_design {};

    iir_prototype butterworth(int num_poles);
    iir_prototype butterworth_shelf(int num_poles, double gain_db);
    iir_prototype chebyshev1(int num_poles, double passband_ripple_db);
    iir_prototype chebyshev2(int num_poles, double stopband_ripple_db);


    struct fir_window
    {
        fir_window(int nfut, int npast);
        const int index0;        ///< Index of the value correspondig to t=0
        Eigen::ArrayXd coeffs;   ///< The filter coefficients, future -> past
        double sum;
    };

    fir_window hann(int nfut, int npast);
    fir_window hamming(int nfut, int npast);
    fir_window blackman(int nfut, int npast);


    struct fir_design
    {
        fir_design() {}

        template <typename Transform>
        fir_design(const fir_window& w, Transform t);

        Eigen::ArrayXf coeffs;
    };


    class lowpass
    {
    public:
        explicit lowpass(double cutoff) : cutoff_(cutoff) {
            if (cutoff<0.0 || cutoff>1.0)
                throw std::invalid_argument("Cutoff frequency must be in [0,1]");
        }
        void operator()(const iir_prototype& p, iir_design& d) const;
        void operator()(const fir_window& w, fir_design& d) const;

    private:
        double cutoff_;
    };

    class highpass
    {
    public:
        explicit highpass(double cutoff) : cutoff_(cutoff) {
            if (cutoff<0.0 || cutoff>1.0)
                throw std::invalid_argument("Cutoff frequency must be in [0,1]");
        }
        void operator()(const iir_prototype& p, iir_design& d) const;
        void operator()(const fir_window& w, fir_design& d) const;

    private:
        double cutoff_;
    };


    class bandpass
    {
    public:
        bandpass(double lo, double hi) : lo_(lo), hi_(hi) {
            if (lo>=hi || lo<0.0 || hi>1.0)
                throw std::invalid_argument("Frequency range must be ordered and in [0,1]");
        }
        void operator()(const iir_prototype& p, iir_design& d) const;
        void operator()(const fir_window& w, fir_design& d) const;

    private:
        double lo_, hi_;
    };


    class bandstop
    {
    public:
        bandstop(double lo, double hi) : lo_(lo), hi_(hi) {
            if (lo>=hi || lo<0.0 || hi>1.0)
                throw std::invalid_argument("Frequency range must be ordered and in [0,1]");
        }
        void operator()(const iir_prototype& p, iir_design& d) const;
        void operator()(const fir_window& w, fir_design& d) const;

    private:
        double lo_, hi_;
    };


    template <typename Transform>
    iir_design::iir_design(const iir_prototype& analog, Transform xform)
    : num_poles_(0), normal_w_(0), normal_g_(1)
    {
        xform(analog, *this);
    }

    template <typename Transform>
    fir_design::fir_design(const fir_window& window, Transform xform)
    {
        xform(window, *this);
    }

}
