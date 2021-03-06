# An exact diagonalization implementation for quantum-dot cellular automata

## Overview

An exact diagonalization implementation for computer simulations of quantum-dot
cellular automata (QCA), a beyond-CMOS computing paradigm ([Original paper proposing
QCA][Lent1993]).

The simulation is written in C++, but exposed as a Python extension module, and
uses [Eigen][] as a linear algebra library and [Boost][] for the Python
integration and unit testing.

## Installation

### Prebuilt binaries

Prebuilt OS X and Linux binary packages for use with the [Anaconda][] Python distribution (or
another distribution using the [conda][] package manager) are available on
[Binstar][]. To use them, first install Anaconda, and then install the qca package
from Binstar.
```
$ conda install -c https://conda.binstar.org/meznom qca
```
This will pull in all dependencies and install everything into the main Anaconda
environment. One notable dependency is boost-python which is also provided by my
Binstar channel.

Alternatively, you can create a new environment and install the qca package
there.
```bash
# Add my Binstar channel permanently to your conda config (~/.condarc)
$ conda config --add channels http://conda.binstar.org/meznom
# Create a new environment 'qca' with a number of useful packages
$ conda create -n qca qca coma ipython ipython-notebook matplotlib numpy scipy
# Activate the 'qca' environment
$ source activate qca
# Do some work, as an example, start the IPython notebook
$ ipython notebook
# Deactivate the 'qca' environment
$ source deactivate
# To remove the 'qca' environment with all packages, you would run:
$ conda remove -n qca --all
```

### Installing from source

Building from source can be a bit tricky. Particularly, I repeatedly had
problems on OS X where it is not uncommon to have multiple versions of Python
and different compilers installed (e.g. when using [MacPorts][] or similar). On
Linux compiling from source is more straightforward, in my experience.

Building requirements:

* [cmake][], version 2.8 or greater.
* [Eigen][], version 3.2.0 or greater.
* [Boost][], version 1.49 or greater
* Compiler and standard library supporting C++11
    * Recent versions of [GCC][] are fine (I am using 4.8)
    * Same for recent versions of [Clang][] (I am using 3.4)
* [Python][] 2.7

Minimal build example:
```
$ mkdir Release
$ cd Release
$ cmake ..
$ make
```

Build example using a MacPorts g++ compiler, MacPorts Boost library, and the
[Anaconda][] Python distribution:
```
$ mkdir Release
$ cd Release
$ CXX=/opt/local/bin/g++-mp-4.8 \
  cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DPYTHON_LIBRARY=~/anaconda/lib/libpython2.7.dylib \
  -DPYTHON_INCLUDE_DIR=~/anaconda/include/python2.7 \
  ..
$ make
```

Build example explicitly setting all configurable paths and using a custom built
Boost library, a custom Eigen installation, and the Anaconda Python
distribution:
```
$ mkdir Release
$ cd Release
$ CXX=/usr/bin/clang++ \
  cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DPYTHON_LIBRARY=~/anaconda/lib/libpython2.7.dylib \
  -DPYTHON_INCLUDE_DIR=~/anaconda/include/python2.7 \
  -DBOOST_ROOT=~/anaconda/ \
  -DEIGEN_INCLUDE_DIR=~/anaconda/include/eigen3/ \
  ..
$ make
```

The `make` command builds the Python extension `_qca.so` which gets copied to
the Python module directory, i.e. `qca/_qca.so`. To build and run the C++ unit
tests:
```
$ cd Release
$ make check
```
To run the Python unit tests:
```
$ python -m unittest qca.test
```
Finally, to install the Python module system wide:
```
$ ./setup.py install
```
Or install the Python module in the home directory:
```
$ ./setup.py install --user
```

It is very important that the Python extension `_qca.so` is linked against the
correct version of the Python and Boost libraries. Especially on OS X this has caused
problems repeatedly. On OS X we can check which libraries are linked with
```
$ otool -L _qca.so
```
and we can change linked libraries with the `install_name_tool`. (On Linux `ldd`
lists linked libraries.) The Boost Python library needs to be linked against the
same Python library (again, we check with otool) and, obviously, that should be
the Python version we are using. Lastly, Boost and this QCA exact
diagonalization code need to be compiled with the same compiler (e.g. Clang
versus GCC) and linked against the same standard library (libc++ versus
libstdc++). If you get obscure compilation errors, that's likely the problem. If
everything compiles correctly, but importing the Python module or running Python
code fails with obscure errors, then likely something is not linked correctly.

Useful cmake variables:

* `BASIS_NUMBER_OF_BITS`. Number of bits for the storage of the basis. Default
  is 32.
* `CMAKE_BUILD_TYPE`. Can be either Release or Debug. Defaults to Release.
* `EIGEN_INCLUDE_DIR`. Where to find Eigen 3. Should be detected automatically.
* `BOOST_ROOT`. Where to find Boost. Should be detected automatically. You can
  also set `BOOST_INCLUDE_DIR` and `BOOST_LIBRARY` directly.
* `USE_CXX_FLAGS_RELEASE`. Compiler flags for Release build type. Reasonable
  defaults should be set automatically, depending on which compiler is used
  (GCC or Clang).
* `USE_CXX_FLAGS_DEBUG`. Compiler flags for Debug build type.

### Conda recipe

[Conda][] recipes for `qca` and its dependencies are available in a [separate
repository][conda-recipes].

## Usage

Usage is via the `qca` python module. The unit tests in `qca.test` provide some
examples.

## License

The code is distributed under the two-clause BSD license. Have a look at the
`LICENSE` file for details.

---
Burkhard Ritter (<burkhard@ualberta.ca>), May 2014


[Lent1993]: http://dx.doi.org/10.1088/0957-4484/4/1/004
[cmake]: http://www.cmake.org
[Eigen]: http://eigen.tuxfamily.org
[Boost]: http://www.boost.org
[GCC]: http://gcc.gnu.org/
[Clang]: http://clang.llvm.org/
[Binstar]: https://binstar.org/meznom
[MacPorts]: https://www.macports.org/
[Anaconda]: https://store.continuum.io/cshop/anaconda/
[Python]: https://www.python.org/
[Conda]: https://github.com/conda/conda
[conda-recipes]: https://bitbucket.org/meznom/conda-recipes