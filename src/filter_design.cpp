#include <disiple/filter_design.hpp>
#include <disiple/fft.hpp>
#include <disiple/impl/next_pow2.hpp>

#include <Eigen/Dense>
#include <unordered_map>
#include <mutex>

namespace disiple {

    void iir_design::add_single(double pole, double zero)
    {
        assert (!(num_poles_&1)); // single comes last

        pairs_.emplace_back(pole, zero);
        ++num_poles_;
    }

    void iir_design::add_conjugate_pair(const complex_t pole, const complex_t zero)
    {
        assert (!(num_poles_&1)); // single comes last

        pairs_.emplace_back(pole, zero, std::conj(pole), std::conj(zero));
        num_poles_ += 2;
    }

    void iir_design::add(const complex_pair& poles, const complex_pair& zeros)
    {
        assert (!(num_poles_&1)); // single comes last
//        assert (poles.is_matched_pair ());
//        assert (zeros.is_matched_pair ());

        pairs_.emplace_back(poles.first, zeros.first, poles.second, zeros.second);
        num_poles_ += 2;
    }

    // notch
    void iir_design::add_notch(double Wo, double BW)
    {
        typedef std::complex<double> complex_t;

        const double gain = 1.0/(1.0+tan(0.5*BW*pi));
        const double w = Wo*pi;

        complex_t zero(cos(w), sin(w));

        Eigen::Matrix2d m(2, 2);
        m << 2.0*gain*zero.real(), -(2.0*gain-1.0), 1.0, 0.0;

        complex_t pole = m.eigenvalues().coeff(0);

        add_conjugate_pair(pole, zero);
    }

    iir_design notch(double Wo, double BW)
    {
        iir_design result;
        result.add_notch(Wo, BW);
        result.set_normal(0.0, 1.0);
        return result;
    }

    // comb
    iir_design comb(int N, double BW)
    {
        iir_design result;

        typedef std::complex<double> complex_t;

        const double gain = 1.0/(1.0+tan(0.25*N*BW*pi));

        Eigen::MatrixXd m(N, N);
        m.fill(0);
        m(0,N-1) = 2*gain-1;
        m.block(1, 0, N-1, N-1).setIdentity();

        complex_t v= m.eigenvalues().coeff(N-1);

        // find poles and zeros
        const int pairs = N / 2;

        for(int i=1;i<=pairs;i++)
        {
            const double W = i* 2*pi/N;
            complex_t zero = complex_t(cos(W), sin(W));
            result.add_conjugate_pair(v.real()*zero, zero);
        }

        if (N & 1)
            result.add_single(v.real(), 1);

        result.set_normal(pi/N, 1.0);

        return result;
    }

    namespace
    {
        template <typename F, typename G>
        fir_window make_window(int nfut, int npast, F f, G g)
        {
            fir_window w(nfut, npast);

            double omega = pi/double(nfut+1);
            for (int t=0; t<nfut; ++t)
                w.sum += w.coeffs[t] = f(omega*double(t+1));
            w.sum += w.coeffs[nfut] = 1.0;
            omega = pi/double(npast+1);
            for (int t=0; t<npast; ++t)
                w.sum += w.coeffs[nfut+1+t] = g(omega*double(t+1));

            return w;
        }

        float fir_gain(Eigen::ArrayXf const& coeffs)
        {
            using pair_type    = std::pair<uint32_t, fft>;
            using storage_type = std::vector<pair_type>;

            // avoid reallocations between runs
            static std::mutex           mut;
            static Eigen::ArrayXf       x;
            static Eigen::ArrayXcf      z;
            static storage_type         ffts;

            std::lock_guard<std::mutex> locker(mut);

            uint32_t n = next_pow2((uint32_t)coeffs.size());
            if ((uint32_t)x.size() != n) {
                x.resize(n);
                z.resize(n/2+1);
            }

            x.segment(0, coeffs.size()) = coeffs;
            x.segment(coeffs.size(), n - coeffs.size()).setZero();

            // check if an fft of length n already exists
            auto it = std::find_if(ffts.begin(), ffts.end(),
                                   [n] (pair_type const& p) { return p.first == n; });

            // if not, construct it
            if (it == ffts.end()) {
                ffts.emplace_back(n, fft(n));
                it = std::prev(ffts.end());
            }

            // perform fft
            it->second(x, z);

            return z.abs().maxCoeff();
        }

        template <typename F>
        void make_filter(const fir_window& w, fir_design& d, double f0, F f)
        {
            const Eigen::DenseIndex n     = w.coeffs.size();
            const Eigen::DenseIndex nfut  = w.index0;
            const Eigen::DenseIndex npast = n - nfut - 1;

            d.coeffs.resize(n);
            d.coeffs[w.index0] = (float)f0;

            for (int t = 1; t <= std::max(nfut, npast); ++t)
            {
                double x = double(t) * pi;
                double y = f(x);
                int i0 = w.index0 - t; if (i0 >= 0) d.coeffs[i0] = float(y * w.coeffs[i0]);
                int i1 = w.index0 + t; if (i1 < n)  d.coeffs[i1] = float(y * w.coeffs[i1]);
            }
        }

        void gain_to_unity(fir_design& d)
        {
            d.coeffs /= fir_gain(d.coeffs);
        }

        void dc_to_zero(fir_design& d, const fir_window& w)
        {
            d.coeffs -= w.coeffs.cast<float>() * float(d.coeffs.sum()/w.sum);
        }
    }

    fir_window::fir_window(int nfut, int npast)
    : index0(nfut), coeffs(nfut+1+npast), sum(0.0)
    {}

    fir_window hann(int nfut, int npast)
    {
        return make_window(nfut, npast,
                           [] (double x) { return 0.5 - 0.5 * cos(x); },
                           [] (double x) { return 0.5 + 0.5 * cos(x); }
                           );
    }

    fir_window hamming(int nfut, int npast)
    {
        return make_window(nfut, npast,
                           [] (double x) { return 0.54 - 0.46 * cos(x); },
                           [] (double x) { return 0.54 + 0.46 * cos(x); }
                           );
    }

    fir_window blackman(int nfut, int npast)
    {
        return make_window(nfut, npast,
                           [] (double x) { return 0.42 - 0.5 * cos(x) + 0.08 * cos(2*x); },
                           [] (double x) { return 0.42 + 0.5 * cos(x) + 0.08 * cos(2*x); }
                           );
    }

    void lowpass::operator()(const fir_window& w, fir_design& d) const
    {
        make_filter(w, d, cutoff_,
                    [&] (double x) { return sin(cutoff_ * x) / x; });

        gain_to_unity(d);
    }

    void highpass::operator()(const fir_window& w, fir_design& d) const
    {
        make_filter(w, d, 1.0 - cutoff_,
                    [&] (double x) { return -sin(cutoff_ * x) / x; });

        dc_to_zero(d, w);
        gain_to_unity(d);
    }

    void bandstop::operator()(const fir_window& w, fir_design& d) const
    {
        make_filter(w, d, 1.0 + lo_ - hi_,
                    [&] (double x) { return (sin(lo_ * x) - sin(hi_ * x)) / x; });

        gain_to_unity(d);
    }

    void bandpass::operator()(const fir_window& w, fir_design& d) const
    {
        make_filter(w, d, lo_ - hi_,
                    [&] (double x) { return (sin(lo_ * x) - sin(hi_ * x)) / x; });

        dc_to_zero(d, w);
        gain_to_unity(d);
    }

}
