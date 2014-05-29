#!/bin/bash

# Ensure we are not using MacPorts, but the native OS X compilers
export PATH=$PREFIX/bin:/bin:/sbin:/usr/bin:/usr/sbin

mkdir -vp ${PREFIX}/lib;
./bootstrap.sh \
    --prefix="${PREFIX}/" \
    --with-libraries=python \
    --with-python=${PYTHON};
./b2 \
    --layout=tagged \
    toolset=clang \
    cxxflags="-std=c++11 -stdlib=libc++ -mmacosx-version-min=10.8" \
    linkflags="-stdlib=libc++ -mmacosx-version-min=10.8" \
    stage;
cp stage/lib/libboost_python* ${PREFIX}/lib/

exit 0
