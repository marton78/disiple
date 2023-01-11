#include <catch2/catch_all.hpp>
#include <disiple/named_params.hpp>

using namespace disiple;

template <int N> struct Foo { static constexpr int foo = N; };
template <int N> struct Bar { static constexpr int bar = N; };
template <typename T> struct Baz { using baz = T; };
template <typename T> struct Qux { using qux = T; };

template <typename... Args>
struct Example
{
    using P = Parameters<
        List<Args...>,
        RequiredType<Baz>,
        OptionalType<Qux, float>,
        RequiredValue<int, Foo>,
        OptionalValue<int, Bar, -1>
    >;
};

TEST_CASE("Named Parameters", "[named_params]")
{
    using Ex1 = Example<Foo<5>, Baz<double>>;
    STATIC_REQUIRE( Ex1::P::foo == 5  );
    STATIC_REQUIRE( Ex1::P::bar == -1 );
    STATIC_REQUIRE( std::is_same<Ex1::P::baz, double>::value );
    STATIC_REQUIRE( std::is_same<Ex1::P::qux, float>::value );

    using Ex2 = Example<Foo<6>, Bar<7>, Baz<int>, Qux<void>>;
    STATIC_REQUIRE( Ex2::P::foo == 6 );
    STATIC_REQUIRE( Ex2::P::bar == 7 );
    STATIC_REQUIRE( std::is_same<Ex2::P::baz, int>::value );
    STATIC_REQUIRE( std::is_same<Ex2::P::qux, void>::value );

    using Ex3 = Example<Bar<7>, Foo<6>, Qux<void>, Baz<int>>;
    STATIC_REQUIRE( Ex3::P::foo == 6 );
    STATIC_REQUIRE( Ex3::P::bar == 7 );
    STATIC_REQUIRE( std::is_same<Ex3::P::baz, int>::value );
    STATIC_REQUIRE( std::is_same<Ex3::P::qux, void>::value );
}
