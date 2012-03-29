project(ExactDiagonalization_tests)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

find_path(
	BOOST_INCLUDE_DIR
	boost/test/unit_test.hpp
	PATHS ~/include ~/src/include ~/o/include ~/o/src/include /opt/local/include /sw/include
)
include_directories(
	${BOOST_INCLUDE_DIR}
	../src
)
find_library(
	BOOST_LIBRARY
	boost_unit_test_framework
	PATHS ~/lib ~/src/lib ~/o/lib ~/o/src/lib /opt/local/lib /sw/lib
)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -USTORAGE_TYPE_OF_FERMIONIC_STATE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -USTORAGE_TYPE_OF_FERMIONIC_STATE")

add_executable(utilitiesTest utilitiesTest.cpp)
target_link_libraries(utilitiesTest ${BOOST_LIBRARY})
add_test(utilities utilitiesTest)

add_executable(basisTest basisTest.cpp)
target_link_libraries(basisTest ${BOOST_LIBRARY})
add_test(basis basisTest)

add_executable(systemTest systemTest.cpp)
target_link_libraries(systemTest ${BOOST_LIBRARY})
add_test(system systemTest)

add_executable(eigenHelpersTest eigenHelpersTest.cpp)
target_link_libraries(eigenHelpersTest ${BOOST_LIBRARY})
add_test(eigenHelpers eigenHelpersTest)

add_executable(qcaTest qcaTest.cpp)
target_link_libraries(qcaTest ${BOOST_LIBRARY})
add_test(qca qcaTest)

add_dependencies(buildtest utilitiesTest basisTest systemTest eigenHelpersTest qcaTest)
