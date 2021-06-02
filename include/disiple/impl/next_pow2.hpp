#pragma once

#include <cstdint>
#include <type_traits>

namespace disiple {

#if defined(_WIN32) && (_MSC_VER >= 1300) && defined(_M_IX86)
#  include <intrin.h>
#  define BUILTIN_CLZ msvc_clz

//from https://github.com/llvm-mirror/libcxx/blob/master/include/support/win32/support.h
__forceinline int msvc_clz(uint64_t mask)
{
    unsigned long where;
    // BitScanReverse scans from MSB to LSB for first set bit.
    // Returns 0 if no set bit is found.
#if defined(_WIN64)
    if (_BitScanReverse64(&where, mask))
        return static_cast<int>(63 - where);
#elif defined(_WIN32)
    // Scan the high 32 bits.
    if (_BitScanReverse(&where, static_cast<unsigned long>(mask >> 32)))
        return static_cast<int>(63 -
                                (where + 32)); // Create a bit offset from the MSB.
    // Scan the low 32 bits.
    if (_BitScanReverse(&where, static_cast<unsigned long>(mask)))
        return static_cast<int>(63 - where);
#else
#error "Implementation of __builtin_clzll required"
#endif
    return 64; // Undefined Behavior.
}

__forceinline int msvc_clz(uint32_t mask)
{
    unsigned long where;
    // Search from LSB to MSB for first set bit.
    // Returns zero if no set bit is found.
    if (_BitScanReverse(&where, mask))
        return static_cast<int>(31 - where);
    return 32; // Undefined Behavior.
}

#elif (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
# define BUILTIN_CLZ gcc_clz
inline __attribute__((always_inline)) int gcc_clz(unsigned int x)       { return __builtin_clz(x); }
inline __attribute__((always_inline)) int gcc_clz(unsigned long x)      { return __builtin_clzl(x); }
inline __attribute__((always_inline)) int gcc_clz(unsigned long long x) { return __builtin_clzll(x); }
#endif

template <int Size, bool Signed>
struct next_pow2_impl
{
    static_assert(!Signed, "Type must be uint32_t or uint64_t.");
};

template <>
struct next_pow2_impl<32, false>
{
    static uint32_t run(uint32_t n)
    {
#ifdef BUILTIN_CLZ
        return n <= 2 ? n : uint32_t(1) << (1 + (31 - BUILTIN_CLZ(n - 1)));
#else
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        return n+1;
#endif
    }
};

template <>
struct next_pow2_impl<64, false>
{
    static uint64_t run(uint64_t n)
    {
#ifdef BUILTIN_CLZ
        return n <= 2 ? n : uint64_t(1) << (1 + (63 - BUILTIN_CLZ(n - 1)));
#else
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n |= n >> 32;
        return n+1;
#endif
    }
};

template <typename T>
T next_pow2(T t) {
    return next_pow2_impl<sizeof(T)*8, std::is_signed<T>::value>::run(t);
}

}
