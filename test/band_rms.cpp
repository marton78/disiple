#include <catch2/catch_all.hpp>
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
    ((typename Scalar, int NChan), Scalar, NChan), 
    (float, Dynamic), (double, Dynamic),
    (float, nchan),   (double, nchan)
) {
    using BPFilter    = IIR<Scalar, Stages<4>, Channels<NChan>>;
    using Expectation = MovingAverage<Scalar, Channels<NChan>>;

    IIRDesign bandpass_design(butterworth(4), Bandpass(.1, .2));
    int expectation_design = 10;

    BandRMS<Scalar, BPFilter, Expectation, Channels<NChan>> brms(bandpass_design, expectation_design);

    BPFilter     bp(bandpass_design);
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

        // Calculate results via BandRMS
        brms.apply(data.col(i), y);

        REQUIRE( (z - y).abs().maxCoeff() <= threshold<Scalar>() );
    }
}
