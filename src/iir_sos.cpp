#include <disiple/impl/iir_sos.hpp>

namespace disiple {

namespace internal {

    using namespace Eigen;

    template <typename Scalar>
    void response(Ref<const Array<Scalar, 4, Dynamic>, Aligned, Stride<4, 1>> coeffs,
                  Scalar                                                      scaling,
                  Ref<const Array<Scalar, Dynamic, 1>>                        f,
                  Ref<Array<std::complex<Scalar>, Dynamic, 1>>                z
                  )
    {
        if (z.size() != f.size())
            throw std::runtime_error("Make sure z has the same size as f.");

        typedef std::complex<Scalar> complex_t;

        Array<complex_t, Dynamic, 1> czn1(f.size()), czn2(f.size());

        czn1.real() = (-Scalar(pi) * f).cos();
        czn1.imag() = (-Scalar(pi) * f).sin();
        czn2.real() = (-Scalar(two_pi) * f).cos();
        czn2.imag() = (-Scalar(two_pi) * f).sin();

        z.fill(complex_t(scaling));

        for (int i=0; i<coeffs.cols(); i++)
        {
            auto c = coeffs.col(i);
            z *= (complex_t(1) + c[0] * czn1 + c[1] * czn2)
               / (complex_t(1) - c[2] * czn1 - c[3] * czn2);
        }
    }

    template <typename Scalar>
    void biquad_from_pz_pair(Ref<Array<Scalar, 4, 1>, Aligned> c,
                             complex_t pole1, complex_t zero1,
                             complex_t pole2, complex_t zero2)
    {
        if (zero1.imag() != 0) { // zero2 == std::conj(zero1))
            c[0] = -Scalar(2 * zero1.real());
            c[1] =  Scalar(std::norm(zero1));
        } else {                // zero2.imag() == 0
            c[0] = -Scalar(zero1.real() + zero2.real());
            c[1] =  Scalar(zero1.real() * zero2.real());
        }

        if (pole1.imag() != 0) { // pole2 == std::conj(pole1)
            c[2] =  Scalar(2 * pole1.real());
            c[3] = -Scalar(std::norm(pole1));
        } else {                 // pole2.imag() == 0
            c[2] =  Scalar(pole1.real() + pole2.real());
            c[3] = -Scalar(pole1.real() * pole2.real());
        }
    }

    // explicitly instantiate for float and double

    template
    void biquad_from_pz_pair<float> (Ref<Array4f, Aligned>,
                                     complex_t, complex_t, complex_t, complex_t);
    template
    void biquad_from_pz_pair<double>(Ref<Array4d, Aligned>,
                                     complex_t, complex_t, complex_t, complex_t);

    template
    void response<float>(Ref<const Array<float, 4, Dynamic>, Aligned, Stride<4, 1>>,
                         float, Ref<const Array<float, Dynamic, 1>>,
                         Ref<Array<std::complex<float>, Dynamic, 1>>);

    template
    void response<double>(Ref<const Array<double, 4, Dynamic>, Aligned, Stride<4, 1>>,
                          double, Ref<const Array<double, Dynamic, 1>>,
                          Ref<Array<std::complex<double>, Dynamic, 1>>);

}

}
