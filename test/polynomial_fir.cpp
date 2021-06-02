#include <catch2/catch.hpp>
#include <disiple/fir.hpp>
#include <disiple/polynomial_fir.hpp>
#include <cstdlib>
#include <Eigen/Cholesky>

using namespace Eigen;
using namespace disiple;

static const size_t nchan = 4;
static const size_t ndata = 1024;

namespace {
    template <typename Scalar> Scalar threshold();
    template <> float  threshold<float>()  { return 1e-4f; }
    template <> double threshold<double>() { return 1e-12; }
}

TEMPLATE_TEST_CASE("Polynomial FIR filter", "[fir_poly]", float, double)
{
    using Scalar = TestType;
    
    typedef Array<Scalar, nchan, 1> arrc;

    // data is channels x time
    Array<Scalar, nchan, Dynamic>  raw_data(nchan, ndata);

    // fill data with signal
    for (size_t t=0; t<ndata; ++t)
        raw_data.col(t) = Scalar(t) * arrc(.1,.2,-.3,.5) + (Scalar)std::rand() / (Scalar)RAND_MAX;

    for (size_t wlen = 128; wlen <= 512; wlen += 128)
    {
        Matrix<Scalar, 2, Dynamic> phi(2, wlen);

        // phi = [ones(1,wlen); wlen:1]
        for (size_t t=0; t<wlen; ++t)
            phi.col(t) = Array<Scalar, 2, 1>(Scalar(1), wlen - t);

        // M = phi * phi' \ phi
        Matrix<Scalar, 2, Dynamic> M = (phi*phi.transpose()).ldlt().solve(phi);
        fir<arrc> f = M.row(1).array();

        // Determined by examining the least squares matrix A*A'\A:
        // When fitting a line ax+b to some values X,
        // the slope of the line is a = w'*X
        // where w = linspace(-w0,w0,n) with n = numel(X) and w0 = 6 / (n + n^2);
        // the increment in w is thus d = (2*w0) / (n-1);
        // linspace(-w0,w0,n) = x * (12 / ((n.^2-1)*n)) - (6*n + 6) / ((n.^2-1)*n)

        Eigen::Array<Scalar, 2, 2> P;
        Eigen::Array<Scalar, 2, 4> Q;

        P <<   6, 6,   // x^0
             -12, 0;   // x^1

        Q << 0, -1, 0, 1,   // x^0
             0, -1, 0, 1;   // x^1

        polynomial_fir<arrc, 2, 2, 4> p(P, Q, wlen);

        Array<Scalar, nchan, Dynamic> y(nchan, ndata), z(nchan, ndata);
        for (size_t t=0; t<ndata; ++t)
        {
            f.apply(raw_data.col(t), y.col(t));
            p.apply(raw_data.col(t), z.col(t));
        }

        Scalar maxdev = (y.block(0,wlen,nchan,ndata-wlen) -
                         z.block(0,wlen,nchan,ndata-wlen)).abs().maxCoeff(); 
        
        REQUIRE( maxdev <= threshold<Scalar>() );
    }

}
