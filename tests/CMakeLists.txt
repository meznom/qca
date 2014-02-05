project(ExactDiagonalization_tests)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(
	../src
)

find_package(Boost 1.49.0 COMPONENTS unit_test_framework)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -USTORAGE_TYPE_OF_FERMIONIC_STATE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -USTORAGE_TYPE_OF_FERMIONIC_STATE")

add_executable(basisTest basisTest.cpp)
target_link_libraries(basisTest ${Boost_LIBRARIES})
add_test(basis basisTest)

add_executable(systemTest systemTest.cpp)
target_link_libraries(systemTest ${Boost_LIBRARIES})
add_test(system systemTest)

add_executable(qcaTest qcaTest.cpp)
target_link_libraries(qcaTest ${Boost_LIBRARIES})
add_test(qca qcaTest)

add_executable(qcaIsingTest qcaIsingTest.cpp)
target_link_libraries(qcaIsingTest ${Boost_LIBRARIES})
add_test(qcaIsing qcaIsingTest)

add_dependencies(
    buildtest
    basisTest
    systemTest
    qcaTest
    qcaIsingTest)
