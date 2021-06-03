cmake_minimum_required (VERSION 3.0)
project (disiple VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set(LIB_SOURCES
    src/fft.cpp
    src/filter_design.cpp
    src/iir_prototype.cpp
    src/iir_sos.cpp
    src/iir_transform.cpp
)

set(TEST_SOURCES
    test/main.cpp
    test/iir.cpp
    test/polynomial_fir.cpp
    test/band_rms.cpp
    test/moving_average.cpp
    test/fir.cpp
    test/next_pow2.cpp
    test/running_stats.cpp
    test/delay.cpp
)

add_library (disiple ${LIB_SOURCES})
target_include_directories(disiple PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(disiple ${CONAN_LIBS})

project (disiple_test)
add_executable(disiple_test ${TEST_SOURCES})
target_include_directories(disiple_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(disiple_test disiple)
target_link_libraries(disiple_test ${CONAN_LIBS})
set_target_properties(disiple_test PROPERTIES OUTPUT_NAME test)
