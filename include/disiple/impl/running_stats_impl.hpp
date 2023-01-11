#pragma once

#include <deque>
#include <Eigen/Core>

namespace disiple {

    /// Filter to calculate running minimum.
    /// Algorithm from: David Lemire, “Streaming Maximum-Minimum Filter Using No
    /// More than Three Comparisons per Element”,
    /// Nordic Journal of Computing,  Vol. 13, 2006

    template <typename Scalar>
    struct RunningMinMaxCoeffs
    {
        RunningMinMaxCoeffs() : len_(0) {}

        explicit RunningMinMaxCoeffs(int len)
        : len_(len) {}

        static Scalar scaling() { return Scalar(1); }
        int length() const { return len_; }

        int len_;
    };

    template <typename Scalar, int Channels, template<typename> class Compare>
    struct RunningMinMaxState
    {
        using Coeffs = RunningMinMaxCoeffs<Scalar>;

        RunningMinMaxState()
        {
            if (Channels != Eigen::Dynamic)
                initialize();
        }

        void setup(const Coeffs& coeffs, int nchans)
        {
            if (nchans != bufs_.size()) {
                bufs_.resize(nchans);
                initialize();
            }
        }

        void initialize()
        {
            for (int i=0; i<bufs_.size(); ++i) {
                bufs_[i].clear();
                bufs_[i].emplace_back(std::numeric_limits<size_t>::max(), Scalar(0));
            }
        }

        template <typename X>
        void apply(const Coeffs& coeffs,
                   const Eigen::ArrayBase<X>& xi, DryRun)
        {
            //David Lemire, “Streaming Maximum-Minimum Filter
            //Using No More than Three Comparisons per Element”
            //Nordic Journal of Computing,  vol. 13, 2006

            Compare<Scalar> cmp;

            for (int c=0; c<bufs_.size(); ++c)
            {
                auto& buf = bufs_[c];
                Scalar x = xi[c];

                //pop oldest element
                if (buf.front().first >= (size_t)coeffs.length())
                    buf.pop_front();

                //remove greater elements from back
                while (!buf.empty() && !cmp(buf.back().second, x))
                    buf.pop_back();

                //increase ages of all remaining elements
                for (auto& p : buf) ++p.first;

                buf.emplace_back(1, x);
            }
        }

        template <typename F, typename Dst>
        void process_result(F f, Eigen::ArrayBase<Dst>& dst) const
        {
            for (int i=0; i<bufs_.size(); ++i)
                f(bufs_[i].front().second, dst[i]);
        }

        template <typename X>
        void apply(const Coeffs& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            apply(coeffs, xi, dry_run);
            process_result([] (Scalar src, Scalar& dst) { dst = src; }, xi);
        }

        using buf_t = std::deque<std::pair<size_t, Scalar>>;
        Eigen::Array<buf_t, Channels, 1> bufs_;
        int length_;
    };


    template <typename Scalar, int Channels>
    struct RunningRangeState
    {
        using Coeffs = RunningMinMaxCoeffs<Scalar>;

        void setup(const Coeffs& coeffs, int nchans)
        {
            smin.setup(coeffs, nchans);
            smax.setup(coeffs, nchans);
        }

        void initialize()
        {
            smin.initialize();
            smax.initialize();
        }

        template <typename X>
        void apply(const Coeffs& coeffs,
                   Eigen::ArrayBase<X>& xi)
        {
            smin.apply(coeffs, xi, dry_run);
            smax.apply(coeffs, xi);
            smin.process_result([] (Scalar src, Scalar& dst) { dst -= src; }, xi);
        }

        template <typename X>
        void apply(const Coeffs& coeffs,
                   const Eigen::ArrayBase<X>& xi, DryRun)
        {
            smin.apply(coeffs, xi, dry_run);
            smax.apply(coeffs, xi, dry_run);
        }

        RunningMinMaxState<Scalar, Channels, std::less>    smin;
        RunningMinMaxState<Scalar, Channels, std::greater> smax;
    };

}
