#include <catch2/catch.hpp>
#include <disiple/impl/next_pow2.hpp>

TEST_CASE("Next power of 2 of a uint32_t", "[next_pow2]")
{
    for (uint8_t k=2; k<31; ++k)
    {
        uint32_t n = uint32_t(1) << k;
        REQUIRE(n == disiple::next_pow2(n));
        REQUIRE(n == disiple::next_pow2(n-1));
        REQUIRE(n == disiple::next_pow2(n/2+1));
    }
}

TEST_CASE("Next power of 2 of a uint64_t", "[next_pow2]")
{
    for (uint8_t k=2; k<63; ++k)
    {
        uint64_t n = uint64_t(1) << k;
        REQUIRE(n == disiple::next_pow2(n));
        REQUIRE(n == disiple::next_pow2(n-1));
        REQUIRE(n == disiple::next_pow2(n/2+1));
    }
}
