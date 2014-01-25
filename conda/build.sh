#!/bin/bash

mkdir Release
cd Release
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBASIS_NUMBER_OF_BITS=32 \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_INSTALL_RPATH=$LD_RUN_PATH \
    -DPYTHON_INCLUDE_PATH:PATH=$PREFIX/include/python${PY_VER} \
    -DPYTHON_LIBRARY:FILEPATH=$PREFIX/lib/libpython${PY_VER}.so \
    ..
make
cd ..

$PYTHON setup.py install

exit 0
