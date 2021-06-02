#include <disiple/fft.hpp>

// these typedefs must come before including Accelerate.h,
// because Accelerate redefines std::complex
typedef std::complex<float>  cf;
typedef std::complex<double> cd;

#if USE_ACCELERATE
#    include <Accelerate/Accelerate.h>
#elif USE_PFFFT
#    include <pffft/pffft.h>
#else
#    include <unsupported/Eigen/FFT>
#endif

namespace disiple {

    static unsigned int ilog2(size_t x) {
        if (x == 0 || (x & (x - 1)))
            throw std::runtime_error("Length must be a power of 2");
        unsigned int y = 0;
        while (x >>= 1) ++y;
        return y;
    }

    fft::fft(size_t length)
    : len_(length), loglen_(ilog2(length))
    {
#if USE_ACCELERATE
        impl_.reset(vDSP_create_fftsetup(loglen_, kFFTRadix2));
        work_.resize(len_);
#elif USE_PFFFT
        impl_.reset(pffft_new_setup(boost::numeric_cast<int>(len_), PFFFT_REAL));
        work_.resize(len_);
#else
        auto* p = new Eigen::FFT<float>;
        p->SetFlag(Eigen::FFT<float>::HalfSpectrum);
        impl_.reset(p);
#endif
    }

    void fft::impl_deleter::operator()(void* p) const
    {
#if USE_ACCELERATE
        vDSP_destroy_fftsetup(static_cast<FFTSetup>(p));
#elif USE_PFFFT
        pffft_destroy_setup(static_cast<PFFFT_Setup*>(p));
#else
        delete static_cast<Eigen::FFT<float>*>(p);
#endif
    }

    using namespace Eigen;

    void fft::operator()(Ref<const ArrayXf, Aligned> x,
                         Ref<ArrayXcf, Aligned> y) const
    {
        if (static_cast<size_t>(x.size()) != len_)
            throw std::runtime_error(std::string("fft: unexpected input size: ") +
                                     std::to_string(x.size()));

        if (static_cast<size_t>(y.size()) != len_/2+1)
            throw std::runtime_error(std::string("fft: unexpected output size: ") +
                                     std::to_string(y.size()));

        //perform fft
#if USE_ACCELERATE
        const size_t halflen = len_/2;
        DSPSplitComplex y_spl = { (float*) y.data(), (float*) y.data() + halflen };
        DSPSplitComplex w_spl = { work_.data(), work_.data() + halflen };

        //copy x's contents to w in split complex format
        vDSP_ctoz(reinterpret_cast<DSPComplex const*>(x.data()), 2, &w_spl, 1, halflen);
        //in-place fft(w), use y as temporary working area
        vDSP_fft_zript(static_cast<FFTSetup>(impl_.get()), &w_spl, 1, &y_spl, loglen_, FFT_FORWARD);
        //un-scale spectrum, see vDSP reference
        vDSP_ztoc(&w_spl, 1, reinterpret_cast<DSPComplex*>(y.data()), 2, halflen);
        y *= 0.5f;
        y[len_/2] = cf(y[0].imag(), 0); y[0] = cf(y[0].real(), 0); //unmix DC and Nyquist
#elif USE_PFFFT
        pffft_transform_ordered(static_cast<PFFFT_Setup*>(impl_.get()),
                                x.data(), reinterpret_cast<float*>(y.data()),
                                work_.data(), PFFFT_FORWARD);
        y[len_/2] = cf(y[0].imag(), 0); y[0] = cf(y[0].real(), 0); //unmix DC and Nyquist
#else
        Eigen::FFT<float>* fft = static_cast<Eigen::FFT<float>*>(impl_.get());
        fft->fwd(y.data(), x.data(), len_);
#endif
    }

}
