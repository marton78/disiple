#include <catch2/catch_all.hpp>
#include <disiple/named_params.hpp>

using namespace disiple;

template <int N> struct Foo { static constexpr int foo = N; };
template <int N> struct Bar { static constexpr int bar = N; };

template <typename... Args>
struct Example
{
    using P = Parameters<
        List<Args...>,
        Required<int, Foo>,
        Optional<int, Bar, -1>
    >;
};

TEST_CASE("Named Parameters", "[named_params]")
{
    using Ex1 = Example<Foo<5>>;
    STATIC_REQUIRE( Ex1::P::foo == 5  );
    STATIC_REQUIRE( Ex1::P::bar == -1 );

    using Ex2 = Example<Foo<6>, Bar<7>>;
    STATIC_REQUIRE( Ex2::P::foo == 6 );
    STATIC_REQUIRE( Ex2::P::bar == 7 );

    using Ex3 = Example<Bar<7>, Foo<6>>;
    STATIC_REQUIRE( Ex3::P::foo == 6 );
    STATIC_REQUIRE( Ex3::P::bar == 7 );
}
