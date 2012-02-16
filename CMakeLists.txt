### generic ###
project(ExactDiagonalization)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

### set up paths and flags ###
find_path(
	EIGEN_INCLUDE_DIR 
    NAMES Eigen/Sparse Eigen/Dense
	PATHS ~/include ~/src/include ~/o/include ~/o/src/include /opt/local/include /sw/include
	PATH_SUFFIXES eigen3 eigen
)
include_directories(
	${EIGEN_INCLUDE_DIR}
)

set(BASIS_NUMBER_OF_BITS "16" CACHE STRING "Number of bits used to store basis. Should be 16, 32 or 64. Generally this should be equal to or larger than 8 times the number of QCA cells.")

### check compiler and set flags appropriately ###
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-debug IS_INTEL_COMPILER)
if (IS_INTEL_COMPILER)
	message(STATUS "Intel compiler detected")
	set(USE_CXX_FLAGS_DEBUG "-O0 -Wall -g" CACHE STRING "")
	#set(USE_CXX_FLAGS_RELEASE "-fast -openmp -Wall -DNDEBUG -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET" CACHE STRING "")
	set(USE_CXX_FLAGS_RELEASE "-fast -Wall -DNDEBUG" CACHE STRING "")
else (IS_INTEL_COMPILER)
	message(STATUS "GNU compiler detected")
	set(USE_CXX_FLAGS_DEBUG "-m64 -O0 -Wall -ggdb" CACHE STRING "")
	set(USE_CXX_FLAGS_RELEASE "-m64 -O3 -Wall -fopenmp -march=native -DNDEBUG" CACHE STRING "")
endif (IS_INTEL_COMPILER)

set(CMAKE_CXX_FLAGS_DEBUG "${USE_CXX_FLAGS_DEBUG} -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET -DSTORAGE_TYPE_OF_FERMIONIC_STATE=uint${BASIS_NUMBER_OF_BITS}_t")
set(CMAKE_CXX_FLAGS_RELEASE "${USE_CXX_FLAGS_RELEASE} -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET -DSTORAGE_TYPE_OF_FERMIONIC_STATE=uint${BASIS_NUMBER_OF_BITS}_t")

# make sure that the default is a Release
if (NOT CMAKE_BUILD_TYPE)
    set (
        CMAKE_BUILD_TYPE Release CACHE STRING 
        "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)


### tests ###
enable_testing()
add_custom_target(buildtest)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check buildtest)
add_subdirectory(tests EXCLUDE_FROM_ALL)

### documentation ###
add_subdirectory(doc EXCLUDE_FROM_ALL)

### our executables ###
add_executable(runQca
    runQca.cpp
    version.hpp
)

# get program version from git
add_custom_command(OUTPUT version.hpp
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND bash updateVersion.hpp.sh
)
