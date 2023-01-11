#pragma once

#include <Eigen/Core>
#include <memory>

namespace disiple {

    /// 1d half-spectrum real-to-complex fft
    class FFT
    {
    public:
        explicit FFT(size_t length);
        void operator()(Eigen::Ref<const Eigen::ArrayXf, Eigen::Aligned>x ,
                        Eigen::Ref<Eigen::ArrayXcf, Eigen::Aligned> y) const;

    private:
        struct ImplDeleter { void operator()(void* p) const; };
        std::unique_ptr<void, ImplDeleter>  impl_;
        size_t                              len_;
        unsigned int                        loglen_;
        mutable Eigen::ArrayXf              work_;
    };

}
