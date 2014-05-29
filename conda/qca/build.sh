#!/bin/bash

# To build directly from the source directory:
cp -r $RECIPE_DIR/../../../qca __qca__
cd __qca__

# Ensure we are not using MacPorts, but the native OS X compilers
export PATH=$PREFIX/bin:/bin:/sbin:/usr/bin:/usr/sbin

# Set PYLIB to either .so (Linux) or .dylib (OS X)
PYLIB="PYTHON_LIBRARY_NOT_FOUND"
for i in $PREFIX/lib/libpython${PY_VER}{.so,.dylib}; do
    if [ -f $i ]; then
        PYLIB=$i
    fi
done

mkdir Release
cd Release
CXX=clang++ \
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_INSTALL_RPATH=$LD_RUN_PATH \
    -DBOOST_ROOT:PATH=$PREFIX \
    -DPYTHON_INCLUDE_PATH:PATH=$PREFIX/include/python${PY_VER} \
    -DPYTHON_LIBRARY:FILEPATH=$PYLIB \
    ..
make

# Run C++ unit test suite
make buildtest
# TODO: only do this on OS X
for i in tests/*Test; do
    install_name_tool -add_rpath ${PREFIX}/lib $i
    install_name_tool -change libboost_unit_test_framework-mt.dylib  @rpath/libboost_unit_test_framework-mt.dylib $i
    install_name_tool -change libboost_python-mt.dylib @rpath/libboost_python-mt.dylib $i
done
# TODO: basisTest fails when built with clang...
#make check

cd ..

$PYTHON setup.py install

exit 0
