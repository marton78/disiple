#pragma once

#include <utility>
#include <optional>
#include <deque>

// TODO: axis AxisElements, not element!
// filter base should accept:
// if axis is regular:
//    scalars
//    axis elements + scalars
// if axis is irregular:
//    axis elements + scalars

namespace disiple {

/// A class representing an axis with irregularly placed elements.
/// @tparam T The type of each value along this axis
template <typename T>
struct IrregularAxis
{
  using Value = T;
};

/// A class representing an axis with regularly spaced elements.
/// The nth value along this axis will be at `first + delta * index`
/// @tparam T The type of each value along this axis
/// @tparam Index The type of the indices
template <typename T, typename Index = int>
struct RegularAxis
{
  using Value = T;
  Value first;       ///< The value of element 0
  Value delta;       ///< Distance between consecutive elements
};

/// A class representing an axis with regularly spaced elements.
/// The nth value along this axis will be at `first + delta * index`
/// @tparam T The type of each value along this axis
/// @tparam Index The type of the indices
template <typename T, typename Index = int>
struct UnitAxis
{
  using Value = T;
};


/// Represents one element along an axis.
template <typename Axis>
struct AxisElement;

template <typename T>
struct AxisElement<IrregularAxis<T>>
{
  using Value = T;
  using Axis  = IrregularAxis<T>;
  Value value;
  Value get_value(Axis const&) const { return value; }
};

template <typename T, typename Index>
struct AxisElement<RegularAxis<T, Index>>
{
  using Value = T;
  using Axis  = RegularAxis<T, Index>;
  Index index;
  Value get_value(Axis const& a) const { return a.first + index * a.delta; }
};

template <typename T, typename Index>
struct AxisElement<UnitAxis<T, Index>>
{
  using Value = T;
  using Axis  = UnitAxis<T, Index>;
  Index index;
  Value get_value(Axis const& a) const { return index; }
};

}
