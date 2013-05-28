#! /bin/bash

PROGRAM_VERSION=`git describe --always --dirty | tr -d "\n"`
if [[ "$PROGRAM_VERSION" == "" ]]; then
    PROGRAM_VERSION="unknown"
fi
if [[ -f version.hpp && "$PROGRAM_VERSION" == "unknown" ]]; then
    exit 0;
fi
sed 's/${PROGRAM_VERSION}/'$PROGRAM_VERSION'/g' version.hpp.in > version.hpp.tmp
if ( diff version.hpp version.hpp.tmp > /dev/null 2>&1 ); then 
   rm version.hpp.tmp
else
    mv version.hpp.tmp version.hpp
fi
