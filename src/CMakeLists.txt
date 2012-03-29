### generic ###
project(ExactDiagonalization_src)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

### precompiled headers ###
#
# I experimented with precompiled headers.
# The file src/precompiledHeaders.hpp would contain dependencies via #include
# statements, for example Eigen.
# 
# It works with gcc, but compilation is not really faster.
# It does not seem to work with clang.
#
#STRING(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" _flags_var_name)
#SET (_args ${${_flags_var_name}})
#
#GET_DIRECTORY_PROPERTY(DIRINC INCLUDE_DIRECTORIES)
#foreach (_inc ${DIRINC})
#    LIST(APPEND _args "-I" ${_inc})
#endforeach(_inc ${DIRINC})
#
#SEPARATE_ARGUMENTS(_args)
#
#add_custom_command(OUTPUT precompiledHeaders.hpp.gch
#    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#    COMMAND ${CMAKE_CXX_COMPILER} ${_args} -x c++-header -c precompiledHeaders.hpp -o ${PROJECT_BINARY_DIR}/precompiledHeaders.hpp.pch # or .pch for clang
#    DEPENDS "${PROJECT_SOURCE_DIR}/precompiledHeaders.hpp"
#    IMPLICIT_DEPENDS CXX "${PROJECT_SOURCE_DIR}/precompiledHeaders.hpp"
#)
#
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -include ${PROJECT_BINARY_DIR}/precompiledHeaders.hpp")
#add_executable(runQca
#    precompiledHeaders.hpp.gch
#    runQca.cpp
#    version.hpp
#)

# get program version from git
add_custom_command(OUTPUT version.hpp
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND bash updateVersion.hpp.sh
)

### our executables ###
add_executable(runQca
    runQca.cpp
    version.hpp
)
