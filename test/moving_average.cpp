#include <catch2/catch.hpp>
#include <disiple/moving_average.hpp>

using namespace Eigen;
using namespace disiple;

static const size_t nchan = 4;

namespace {
    template <typename Scalar> Scalar threshold();
    template <> float  threshold<float>()  { return 1e-5f; }
    template <> double threshold<double>() { return 1e-10; }
    template <> int    threshold<int>()    { return 0; }
}

TEMPLATE_TEST_CASE_SIG("Moving Average Filter", "[running_stats]",
    ((typename Scalar, int Dim), Scalar, Dim), 
    (float, Dynamic), (double, Dynamic), (int, Dynamic),
    (float, nchan),   (double, nchan),   (int, nchan) 
) {
    const int W = 10; // Window length of the filter

    disiple::moving_average<Array<Scalar, Dim, 1>, Dynamic, 1> rmean1(W);
    disiple::moving_average<Array<Scalar, Dim, 1>, Dynamic, 2> rmean2(W);
    disiple::moving_average<Array<Scalar, Dim, 1>, Dynamic, 3> rmean3(W);

    Array<Scalar, nchan, Dynamic> data = (ArrayXXf::Random(nchan, 97) * 10 + 20).cast<Scalar>();
    Array<Scalar, nchan, 1> y1, y2, y3, z;

    for (int i=0; i<data.cols(); ++i)
    {
        // Calculate results manually
        auto block = i<W ? data.block(0,     0, nchan, i+1)
                         : data.block(0, i-W+1, nchan,   W);
        z = block.rowwise().mean();

        rmean1.apply(data.col(i), y1);
        
        // Check that mavg<1>(data) average equals manual result
        REQUIRE( (z - y1).abs().maxCoeff() <= threshold<Scalar>() );

        rmean2.apply(y1, y2);
        rmean3.apply(data.col(i), y3);

        // Check that mavg<2>(mavg<1>(data)) equals mavg<3>(data)
        REQUIRE( (y2 - y3).abs().maxCoeff() <= threshold<Scalar>() );
    }
}
