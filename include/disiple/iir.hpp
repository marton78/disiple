#pragma once

#include <disiple/impl/filter_base.hpp>
#include <disiple/impl/iir_sos.hpp>
#include <disiple/impl/iir_impl.hpp>
#include <disiple/filter_design.hpp>
#include <disiple/named_params.hpp>

namespace disiple {

    template <IIRImplementation N>
    struct Implementation { static constexpr IIRImplementation implementation = N; };

    /// Infinite impulse response (IIR) digital filter
    template <typename Scalar, typename... Options>
    class IIR : public FilterBase<Scalar, IIR<Scalar, Options...>>,
                public Parameters<
                            List<Options...>,
                            OptionalValue<IIRImplementation, Implementation, DF2T>,
                            OptionalValue<int, Stages, Eigen::Dynamic>,
                            OptionalValue<int, Channels, 1>
                        >
    {
    public:
        using State = IIRImpl<Scalar, IIR::stages, IIR::channels, IIR::implementation>;
        using Coeffs = SecondOrderSections<Scalar, IIR::stages>;

        IIR() {}
        IIR(IIRDesign const& d) : coeffs_(d) {}

        int num_stages() const { return coeffs_.num_stages(); }
        
    private:
        friend FilterBase<Scalar, IIR<Scalar, Options...>>;
        State  state_;
        Coeffs coeffs_;
    };

}
