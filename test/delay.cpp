#include <catch2/catch.hpp>
#include <disiple/delay.hpp>
#include <type_traits>

using namespace Eigen;
using namespace disiple;

static const size_t nchan = 1;
static const size_t ndata = 25;

template <typename Scalar>
struct TestFixtureDelay
{
    enum { Length = 5 };

    // data is channels x time
    Array<Scalar, Dynamic, Dynamic>  raw_data, ref_data;

    TestFixtureDelay()
    : raw_data(nchan, ndata)
    , ref_data(nchan, ndata)
    {
        // fill data with sawtooth signal
        for (size_t t=0; t<ndata; ++t)
            for (size_t c=0; c<nchan; ++c)
                raw_data(c,t) = (t + c) % 16;

        ref_data.leftCols(Length).setConstant(Scalar(0));
        ref_data.rightCols(ndata-Length) = raw_data.leftCols(ndata-Length);
    }
};

TEMPLATE_TEST_CASE_SIG("Delay filter", "[delay]",
    ((typename Scalar, int NChan, bool DynLength), Scalar, NChan, DynLength), 
    (float,  Dynamic, true),  (double, Dynamic, true),  (int, Dynamic, true), 
    (float,  nchan,   true),  (double, nchan,   true),  (int, nchan,   true), 
    (float,  Dynamic, false), (double, Dynamic, false), (int, Dynamic, false), 
    (float,  nchan,   false), (double, nchan,   false), (int, nchan,   false)
) {
    TestFixtureDelay<Scalar> fix;
    enum { Length = TestFixtureDelay<Scalar>::Length };
    Array<Scalar, Dynamic, Dynamic> y(nchan, ndata);
    using Filter = delay<Array<Scalar, NChan, 1>, DynLength ? Dynamic : Length>;

    SECTION("in_place", "Delay filter applied in place") {
        Filter f(Length);
        y = fix.raw_data;
        f.apply(y.block(0,  0, nchan,       10));
        f.apply(y.block(0, 10, nchan, ndata-10));
        Scalar maxdev = (fix.ref_data - y).abs().maxCoeff();
        REQUIRE( maxdev == 0 );
    }

    SECTION("to_other", "FIR filter applied into another array") {
        Filter f(Length);
        f.apply(fix.raw_data.block(0,  0, nchan,       10), y.block(0,  0, nchan,       10));
        f.apply(fix.raw_data.block(0, 10, nchan, ndata-10), y.block(0, 10, nchan, ndata-10));
        Scalar maxdev = (fix.ref_data - y).abs().maxCoeff();
        REQUIRE( maxdev == 0 );
    }

    SECTION("update_only", "FIR filter state updated without applying it") {
        Filter f(Length);
        y = fix.raw_data;
        f.apply(y.block(0,  0, nchan,       10), dry_run);
        f.apply(y.block(0, 10, nchan, ndata-10));
        Scalar maxdev0 = (fix.raw_data.block(0,  0, nchan,       10) - y.block(0,  0, nchan,       10)).abs().maxCoeff();
        Scalar maxdev1 = (fix.ref_data.block(0, 10, nchan, ndata-10) - y.block(0, 10, nchan, ndata-10)).abs().maxCoeff();
        REQUIRE( maxdev0 == 0 );
        REQUIRE( maxdev1 == 0 );
    }

}

