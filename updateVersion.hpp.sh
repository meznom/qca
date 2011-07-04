#! /bin/bash

GIT_PROGRAM_VERSION=`git describe --always | tr -d "\n"`
sed 's/${GIT_PROGRAM_VERSION}/'$GIT_PROGRAM_VERSION'/g' version.hpp.in > version.hpp.tmp
if ( diff version.hpp version.hpp.tmp > /dev/null 2>&1 ); then 
   rm version.hpp.tmp
else
    mv version.hpp.tmp version.hpp
fi
