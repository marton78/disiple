# disiple - Digital Signal Processing Library for Eigen

I wrote this code during my time at Nielsen Telemedical. It is hereby released to the public
with the kind permission of [TeleMedi GmbH](https://www.telemedi.tech/).

The code was originally based on Vinnie Falco's [DSPFilters](https://github.com/vinniefalco/DSPFilters)
library, but has been heavily refactored and extended. Therefore, it still builds "on the
work of cherished luminaries such as Sophocles Orfanidis, Andreas Antoniou, Martin Holters,
and Udo Zolzer" (for the IIR filter design functions), but also on [Daniel Lemire](https://lemire.me/)'s ingenious
paper "Streaming Maximum-Minimum Filter Using No More than Three Comparisons per Element".
My work on [FIR filtering with polynomial coefficients](docs/polynomial_filtering.pdf) is
-- to the best of my knowledge -- novel, but it might well be a rediscovery of something ancient.

The library contains:

* IIR filtering based on second order sections, implemented as:
    * Direct Form 1
    * Direct Form 2
    * Direct Form 2 Transposed
* FIR filtering
* A moving average filter, optionally with multiple stages
* An efficient, recursive implementation of FIR filtering with polynomial coefficients,
  useful e.g. for peak detection by fitting a polynomial to a time series
* Daniel Lemire's efficient streaming Maximum-Minimum Filter
* IIR filter design prototypes: Butterworth, Butterworth shelf, Chebyshev Type 1 and 2
* FIR filter design prototypes: Hann, Hamming and Blackman
* Transformations of prototypes to Lowpass, Highpass, Bandpass and Bandstop IIR and FIR filters
* All code is implemented using the [Eigen](https://eigen.tuxfamily.org/) C++ template library for linear algebra, benefitting from its SIMD operations.

## Building

To build the tests, you'll need the [Conan](https://conan.io/) package manager and [CMake](https://cmake.org/).

```
mkdir build
cd build
conan install ..
cmake
cmake --build .
./bin/test
```