#include <catch2/catch_all.hpp>
#include <disiple/running_statistics.hpp>

using namespace Eigen;
using namespace disiple;

static const size_t nchan = 4;

namespace {

    template <typename Scalar> Scalar threshold();
    template <> float  threshold<float>()  { return 1e-5f; }
    template <> double threshold<double>() { return 1e-10; }
    template <> int    threshold<int>()    { return 0; }

}

TEMPLATE_TEST_CASE_SIG("Running Min, Max and Range", "[running_stats]",
    ((typename Scalar, int NChan), Scalar, NChan), 
    (float, Dynamic), (double, Dynamic), (int, Dynamic),
    (float, nchan),   (double, nchan),   (int, nchan) 
) {
    const int W = 10; // Window length of the filter

    disiple::RunningMin  <Scalar, Channels<NChan>> fmin(W);
    // disiple::RunningMax  <Scalar, Channels<NChan>> fmax(W);
    // disiple::RunningRange<Scalar, Channels<NChan>> frng(W);

    Array<Scalar, nchan, Dynamic> data = (ArrayXXf::Random(nchan, 97) * 10 + 20).cast<Scalar>();
    Array<Scalar, nchan, 1> ymin, ymax, yrng, zmin, zmax, zrng;

    for (int i=0; i<data.cols(); ++i)
    {
        // Calculate results manually
        auto block = i<W ? data.block(0,     0, nchan, i+1)
                         : data.block(0, i-W+1, nchan,   W);
        zmin = block.rowwise().minCoeff();
        zmax = block.rowwise().maxCoeff();
        zrng = zmax - zmin;

        // Calculate results using the filters
        fmin.apply(data.col(i), ymin);
        // fmax.apply(data.col(i), ymax);
        // frng.apply(data.col(i), yrng);

        REQUIRE( (zmin - ymin).abs().maxCoeff() <= threshold<Scalar>() );
        // REQUIRE( (zmax - ymax).abs().maxCoeff() <= threshold<Scalar>() );
        // REQUIRE( (zrng - yrng).abs().maxCoeff() <= threshold<Scalar>() );
    }
}
