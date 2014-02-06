# QCA exact diagonalization

## Installation

Requirements:

* [cmake][1], version 2.8 or greater.
* [Eigen][2], version 3.2.0 or greater.
* [Boost][3], version 1.49 or greater (1.48 does not seem to work)
* works with gcc (I am using 4.8) and clang (I am using 3.0)
* Python 2.7 for the Python bindings

Build, for example:

    $ mkdir Release
    $ cd Release
    $ rm -f CMakeCache.txt && \
      CXX=/opt/local/bin/g++-mp-4.8 \
      cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DPYTHON_LIBRARY=~/anaconda/lib/libpython2.7.dylib \
      -DPYTHON_INCLUDE_DIR=~/anaconda/include/python2.7 \
      ..
    # Builds the _qca.so C++-Python extension
    $ make
    # C++ tests
    $ make check
    $ cd ..
    # Python tests
    $ python -m unittest qca.test
    $ ./setup.py install --user

This clears the cmake cache and builds the program with the GCC compiler in
Release mode and with a 32 bits basis (the default).

It is very important that the Python extension `_qca.so` is linked against the
proper version of the Python library. Especially on OS X this has caused
problems repeatedly. On OS X we can check which libraries are linked with

    $ otool -L _qca.so

and we can change linked libraries with the `install_name_tool`. (On Linux `ldd`
lists linked libraries.) The Boost Python library needs to be linked against the
same Python library (again, we check with otool) and, obviously, that should be
the Python version we are using. Lastly, I've had issues when the Boost Python
library and this QCA exact diagonalization code were not compiled with the same
compiler (e.g. Clang versus GCC).

Useful cmake variables:

* `BASIS_NUMBER_OF_BITS`. Number of bits for the storage of the basis. Default
  is 32.
* `CMAKE_BUILD_TYPE`. Can be either Release or Debug. Defaults to Release.
* `EIGEN_INCLUDE_DIR`. Where to find Eigen 3. Should be detected automatically.
* `BOOST_INCLUDE_DIR`. Where to find Boost. Should be detected automatically.
* `BOOST_LIBRARY`. Where to find the boost dynamic library. Should be detected
  automatically.
* `USE_CXX_FLAGS_RELEASE`. Compiler flags for Release build type. Reasonable
  defaults should be set automatically, depending on which compiler is used
  (GCC or Clang).
* `USE_CXX_FLAGS_DEBUG`. Compiler flags for Debug build type.

## Building a Conda package

The QCA exact diagonalization module now comes with a [conda][4] recipe, so
conda can be used to build and install a binary package. In the `source_dir`:

    $ conda build ./conda/
    $ conda install --use-local qca

Note that the conda recipe fetches the latest version of this QCA exact
diagonalization module from Bitbucket / git.

## Usage

Usage is via the `qca` python module. The unit tests in `qca.test` provide some
examples.

---
Burkhard Ritter (<burkhard@ualberta.ca>), February 2014


[1]: http://www.cmake.org
[2]: http://eigen.tuxfamily.org
[3]: http://www.boost.org
[4]: https://store.continuum.io/cshop/anaconda/
