#include <disiple/filter_design.hpp>

#if defined(_MSC_VER) && _MSC_VER < 1800
inline double asinh(double x) { return std::log(x + std::sqrt(x*x + 1.0)); }
#endif

namespace disiple {

    iir_prototype butterworth(int num_poles)
    {
        iir_prototype result;
        result.set_normal(0, 1);

        const double n2 = 2 * num_poles;
        const int pairs = num_poles / 2;

        for (int i = 0; i < pairs; ++i)
        {
            double theta = pi_half + (2 * i + 1) * pi/n2;
            complex_t c(cos(theta), sin(theta));
            result.add_conjugate_pair(c, std::numeric_limits<double>::infinity());
        }

        if (num_poles & 1)
            result.add_single(-1, std::numeric_limits<double>::infinity());

        return result;
    };

    iir_prototype butterworth_shelf(int num_poles, double gain_db)
    {
        iir_prototype result;
        result.set_normal(pi, 1);

        const double n2 = 2 * num_poles;
        const double gz = -pow(10., gain_db/20/n2);
        const double gp = 1.0 / gz;
        const int pairs = num_poles / 2;

        for (int i = 1; i <= pairs; ++i)
        {
            double theta = pi_half - (2 * i - 1) * pi/n2;
            complex_t c(cos(theta), sin(theta));
            result.add_conjugate_pair(gp * c, gz * c);
        }

        if (num_poles & 1)
            result.add_single(gp, gz);

        return result;
    }


    // Chebyshev Type I
    // http://cnx.org/content/m16906/latest/

    iir_prototype chebyshev1(int num_poles, double ripple_db)
    {
        iir_prototype result;

        const double e2 = pow(10.0, ripple_db * 0.1);
        const double eps = std::sqrt(e2 - 1.);
        const double v0 = asinh(1. / eps) / num_poles;
        const double sinh_v0 = -sinh(v0);
        const double cosh_v0 = cosh(v0);
        const double n2 = 2 * num_poles;
        const int pairs = num_poles / 2;

        for (int i = 0; i < pairs; ++i)
        {
            double theta = (2*i+1-num_poles) * pi / n2;
            double a = sinh_v0 * cos(theta);
            double b = cosh_v0 * sin(theta);

            result.add_conjugate_pair(complex_t(a, b), std::numeric_limits<double>::infinity());
        }

        if (num_poles & 1) {
            result.add_single(sinh_v0, std::numeric_limits<double>::infinity());
            result.set_normal(0, 1);
        } else
            result.set_normal(0, pow(10, -ripple_db/20.));

        return result;
    }


    // Chebyshev Type II

    iir_prototype chebyshev2(int num_poles, double ripple_db)
    {
        iir_prototype result;
        result.set_normal(0.0, 1.0);

        const double eps = std::sqrt((pow(10.0, ripple_db * 0.1)-1));
        const double v0 = asinh(eps) / num_poles;
        const double sinh_v0 = -sinh(v0);
        const double cosh_v0 = cosh(v0);
        const double fn = pi_half/num_poles;

        int k = 1;
        for (int i = num_poles / 2; --i >= 0; k+=2)
        {
            double theta = (k - num_poles) * fn;
            double a = sinh_v0 * cos(theta);
            double b = cosh_v0 * sin(theta);
            double d2 = a*a + b*b;
            double im = 1.0/cos(k*fn);
            result.add_conjugate_pair(complex_t(a/d2, b/d2), complex_t(0, im));
        }

        if (num_poles & 1)
            result.add_single(1.0/sinh_v0, std::numeric_limits<double>::infinity());

        return result;
    }

}
