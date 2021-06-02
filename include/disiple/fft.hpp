#pragma once

#include <Eigen/Core>
#include <memory>

namespace disiple {

    /// 1d half-spectrum real-to-complex fft
    class fft
    {
    public:
        explicit fft(size_t length);
        void operator()(Eigen::Ref<const Eigen::ArrayXf, Eigen::Aligned>x ,
                        Eigen::Ref<Eigen::ArrayXcf, Eigen::Aligned> y) const;

    private:
        struct impl_deleter { void operator()(void* p) const; };
        std::unique_ptr<void, impl_deleter> impl_;
        size_t                              len_;
        unsigned int                        loglen_;
        mutable Eigen::ArrayXf              work_;
    };

}
