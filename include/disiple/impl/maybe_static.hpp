#pragma once

#include <Eigen/Core>

namespace disiple {

    template <int N, typename Tag = void>
    struct MaybeStatic {
        enum { StaticValue = N };
        MaybeStatic() = default;
        MaybeStatic(int value) { assert(value == StaticValue); }
        constexpr bool is_dynamic() const noexcept { return false; }
        constexpr int get() const noexcept { return StaticValue; }
    };

    template <typename Tag>
    struct MaybeStatic<Eigen::Dynamic, Tag> {
        enum { StaticValue = Eigen::Dynamic };
        MaybeStatic(int value) : dynamic_(value) {}
        constexpr bool is_dynamic() const noexcept { return true; }
        int get() const noexcept { return dynamic_; }
    private:
        int dynamic_;
    };

}
