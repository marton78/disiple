#include <catch2/catch.hpp>
#include <disiple/band_rms.hpp>
#include <disiple/iir.hpp>
#include <disiple/moving_average.hpp>

using namespace Eigen;
using namespace disiple;

static const size_t nchan = 4;

namespace {

    template <typename Scalar> Scalar threshold();
    template <> float  threshold<float>()  { return 1e-5f; }
    template <> double threshold<double>() { return 1e-10; }

}

TEMPLATE_TEST_CASE_SIG("Running Bandpass Root Mean Square", "[running_stats]",
    ((typename Scalar, int Dim), Scalar, Dim), 
    (float, Dynamic), (double, Dynamic),
    (float, nchan),   (double, nchan)
) {
    using Element     = Array<Scalar, Dim, 1>;
    using Bandpass    = iir<Element, 4>;
    using Expectation = moving_average<Element>;

    iir_design bandpass_design(butterworth(4), bandpass(.1, .2));
    int expectation_design = 10;

    band_rms<Element, Bandpass, Expectation> brms(bandpass_design, expectation_design);

    Bandpass     bp(bandpass_design);
    Expectation  ex(expectation_design);

    Array<Scalar, nchan, Dynamic> data = (ArrayXXf::Random(nchan, 97) * 10 + 20).cast<Scalar>();
    Array<Scalar, nchan, 1> z, y;

    for (int i=0; i<data.cols(); ++i)
    {
        // Calculate results manually
        bp.apply(data.col(i), z);     // bandpass
        z *= z;                       // square
        ex.apply(z);                  // mean
        z = (z * Scalar(2)).sqrt();   // root

        // Calculate results via band_rms
        brms.apply(data.col(i), y);

        REQUIRE( (z - y).abs().maxCoeff() <= threshold<Scalar>() );
    }
}
