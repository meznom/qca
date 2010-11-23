### generic ###
project(diagonaliser)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)


### set up paths and flags ###
set(EIGEN_INCLUDE_DIR "~/src/include" CACHE PATH "Include path for eigen2")
include_directories(
	${EIGEN_INCLUDE_DIR}
)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -ggdb")
# TODO: why is there no -march=native on mac?
# TODO: use newer gcc on mac
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fopenmp -mtune=native -DNDEBUG")

# make sure that the default is a Release
if (NOT CMAKE_BUILD_TYPE)
    set (
        CMAKE_BUILD_TYPE Release CACHE STRING 
        "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)

### our executables ###
add_executable(runDiagonaliser
    runDiagonaliser.cpp 
)

add_executable(runQca
    runQca.cpp
)
