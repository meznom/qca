# QCA exact diagonalization

## Installation

Requirements:

* [cmake][1], version 2.8 or greater.
* [Eigen][2], version 3.1.3 or greater.
* [Boost][3], version 1.49 or greater (1.48 does not seem to work)
* [Doxygen][4] to generate the documentation.
* works with gcc (I am using 4.8) and clang (I am using 3.0)
* Python 2.7 for the Python bindings

Build:

    $ mkdir build_dir
    $ cd build_dir
    $ ccmake source_dir
    $ make
    # To build and run the unit tests:
    $ make check
    # To only build the tests:
    $ make buildtest
    # To generate the documentation:
    $ make doc
    # To install the Python module
    $ cd ..
    $ ./setup.py install --user

For example, in the `source_dir`:

    $ mkdir Release
    $ cd Release
    $ rm -f CMakeCache.txt && CXX=/opt/local/bin/clang++-mp-3.0 cmake \
      -DCMAKE_BUILD_TYPE=Release -DBASIS_NUMBER_OF_BITS=32 \
      -DPYTHON_LIBRARY="/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Python" ..
    $ make
    $ cd ..
    $ ./setup.py install --user

This clears the cmake cache and builds the program with the Clang compiler in
Release mode and with a 32 bits basis.

Useful cmake variables:

* `BASIS_NUMBER_OF_BITS`. Number of bits for the storage of the basis. Default
  is 16, but for larger systems (e.g. 3 cells with the grand canonical model)
  32 bits is required.
* `CMAKE_BUILD_TYPE`. Can be either Release or Debug. Defaults to Release.
* `EIGEN_INCLUDE_DIR`. Where to find Eigen 3. Should be detected automatically.
* `BOOST_INCLUDE_DIR`. Where to find Boost. Should be detected automatically.
* `BOOST_LIBRARY`. Where to find the boost dynamic library. Should be detected
  automatically.
* `USE_CXX_FLAGS_RELEASE`. Compiler flags for Release build type. Reasonable
  defaults should be set automatically, depending on which compiler is used
  (GCC or Clang).
* `USE_CXX_FLAGS_DEBUG`. Compiler flags for Debug build type.

## Usage

TODO: This section is out of date.

After compilation a file runQca should be in the `build_dir`. Run `./src/runQca
--help` to get started.

For example:

    $ ./src/runQca 'model: grandcanonical, t: 1, td: 0, V0: 1000, \
      layout: {type: wire, epc: 2, a: 0.01, b: 1000}, observables: {E: yes}'

Calculates the energy spectrum for the one-cell-wire grand-canonical system,
with 2 electrons per cell (epc) and the driving cell placed very, very far away
(b=1000).

    $ ./src/runQca 'model: fixed, t: 1, td: 0, V0: 1000, beta: [0.01,0.1,1,10], \
      layout: {type: wire, epc: 2, V1: 100, boa: [1.75,2.5], Pext: 1}, \
      changing: [layout.boa, beta], observables: {P: all}'

Computes the polarization of all cells over beta, for the two cell fixed charge
system and for two different values of boa (b over a, that is, b in units of
a). b is the inter-cell spacing.

## Performance

I am testing the performance of clang 3.0, gcc 4.6.3, and gcc 4.7.0. 

    $ /opt/local/bin/clang++-mp-3.0 --version
    clang version 3.0 (tags/RELEASE_30/final)
    $ /opt/local/bin/g++-mp-4.6 --version
    g++-mp-4.6 (MacPorts gcc46 4.6.3_2) 4.6.3
    $ /opt/local/bin/g++-mp-4.7 --version
    g++-mp-4.7 (MacPorts gcc47 4.7.0_3) 4.7.0
    $ ( time nice -19 ./Release-clang/src/runQca \
            'model: fixed, V0: 1000, beta: 1, t: 1, layout: {cells: 3, a: 0.01, \
             b: 0.02, Pext: 1, epc: 2}, observables: {P: all, P2: all, N: all}' \
        && \
        time nice -19 ./Release-gcc/src/runQca \
            'model: fixed, V0: 1000, beta: 1, t: 1, layout: {cells: 3, a: 0.01, \
             b: 0.02, Pext: 1, epc: 2}, observables: {P: all, P2: all, N: all}' \
        && \
        time nice -19 ./Release-gcc47/src/runQca \
            'model: fixed, V0: 1000, beta: 1, t: 1, layout: {cells: 3, a: 0.01, \
             b: 0.02, Pext: 1, epc: 2}, observables: {P: all, P2: all, N: all}' \
      ) > output 2>&1

From this I get:

    clang 3.0: 25m37.286s
    gcc 4.6.3: 25m33.630s
    gcc 4.7.0: 30m6.368s

So clang 3.0 and gcc 4.6.3 are essentially on par while gcc 4.7.0 is
considerably slower.

The compilation time with clang is a little bit more than two times faster than
with gcc 4.6.3 (~10 seconds instead of ~25 seconds).


Burkhard Ritter (<burkhard@ualberta.ca>), May 2013


[1]: http://www.cmake.org
[2]: http://eigen.tuxfamily.org
[3]: http://www.boost.org
[4]: http://www.doxygen.org
