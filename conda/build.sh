#!/bin/bash

# Set PYLIB to either .so (Linux) or .dylib (OS X)
PYLIB="PYTHON_LIBRARY_NOT_FOUND"
for i in $PREFIX/lib/libpython${PY_VER}{.so,.dylib}; do
    if [ -f $i ]; then
        PYLIB=$i
    fi
done

mkdir Release
cd Release
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_INSTALL_RPATH=$LD_RUN_PATH \
    -DPYTHON_INCLUDE_PATH:PATH=$PREFIX/include/python${PY_VER} \
    -DPYTHON_LIBRARY:FILEPATH=$PYLIB \
    ..
make
cd ..

$PYTHON setup.py install

exit 0
