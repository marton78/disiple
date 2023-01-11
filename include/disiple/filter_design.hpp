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

    struct IIRPrototype;

    //
    // Describes a digital filter as a collection of poles and zeros along with
    // normalization information to achieve a specified gain at a specified
    // frequency. The poles and zeros lie in the z plane.
    //

    class IIRDesign
    {
    public:
        IIRDesign() : num_poles_(0) { pairs_.reserve(16); }

        template <typename Transform>
        IIRDesign(const IIRPrototype& p, Transform t);

        void add_single(double pole, double zero);
        void add_conjugate_pair(const Complex pole, const Complex zero);
        void add(const ComplexPair& poles, const ComplexPair& zeros);
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

    IIRDesign notch(double Wo, double BW);
    IIRDesign comb(int N, double BW);

    struct IIRPrototype : IIRDesign {};

    IIRPrototype butterworth(int num_poles);
    IIRPrototype butterworth_shelf(int num_poles, double gain_db);
    IIRPrototype chebyshev1(int num_poles, double passband_ripple_db);
    IIRPrototype chebyshev2(int num_poles, double stopband_ripple_db);


    struct FIRWindow
    {
        FIRWindow(int nfut, int npast);
        const int index0;        ///< Index of the value correspondig to t=0
        Eigen::ArrayXd coeffs;   ///< The filter coefficients, future -> past
        double sum;
    };

    FIRWindow hann(int nfut, int npast);
    FIRWindow hamming(int nfut, int npast);
    FIRWindow blackman(int nfut, int npast);


    struct FIRDesign
    {
        FIRDesign() {}

        template <typename Transform>
        FIRDesign(const FIRWindow& w, Transform t);

        Eigen::ArrayXf coeffs;
    };


    class Lowpass
    {
    public:
        explicit Lowpass(double cutoff) : cutoff_(cutoff) {
            if (cutoff<0.0 || cutoff>1.0)
                throw std::invalid_argument("Cutoff frequency must be in [0,1]");
        }
        void operator()(const IIRPrototype& p, IIRDesign& d) const;
        void operator()(const FIRWindow& w, FIRDesign& d) const;

    private:
        double cutoff_;
    };

    class Highpass
    {
    public:
        explicit Highpass(double cutoff) : cutoff_(cutoff) {
            if (cutoff<0.0 || cutoff>1.0)
                throw std::invalid_argument("Cutoff frequency must be in [0,1]");
        }
        void operator()(const IIRPrototype& p, IIRDesign& d) const;
        void operator()(const FIRWindow& w, FIRDesign& d) const;

    private:
        double cutoff_;
    };


    class Bandpass
    {
    public:
        Bandpass(double lo, double hi) : lo_(lo), hi_(hi) {
            if (lo>=hi || lo<0.0 || hi>1.0)
                throw std::invalid_argument("Frequency range must be ordered and in [0,1]");
        }
        void operator()(const IIRPrototype& p, IIRDesign& d) const;
        void operator()(const FIRWindow& w, FIRDesign& d) const;

    private:
        double lo_, hi_;
    };


    class Bandstop
    {
    public:
        Bandstop(double lo, double hi) : lo_(lo), hi_(hi) {
            if (lo>=hi || lo<0.0 || hi>1.0)
                throw std::invalid_argument("Frequency range must be ordered and in [0,1]");
        }
        void operator()(const IIRPrototype& p, IIRDesign& d) const;
        void operator()(const FIRWindow& w, FIRDesign& d) const;

    private:
        double lo_, hi_;
    };


    template <typename Transform>
    IIRDesign::IIRDesign(const IIRPrototype& analog, Transform xform)
    : num_poles_(0), normal_w_(0), normal_g_(1)
    {
        xform(analog, *this);
    }

    template <typename Transform>
    FIRDesign::FIRDesign(const FIRWindow& window, Transform xform)
    {
        xform(window, *this);
    }

}
