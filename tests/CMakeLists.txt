project(ExactDiagonalization_tests)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(
	../src
)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -USTORAGE_TYPE_OF_FERMIONIC_STATE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -USTORAGE_TYPE_OF_FERMIONIC_STATE")

add_executable(utilitiesTest utilitiesTest.cpp)
target_link_libraries(utilitiesTest ${Boost_LIBRARIES})
add_test(utilities utilitiesTest)

add_executable(basisTest basisTest.cpp)
target_link_libraries(basisTest ${Boost_LIBRARIES})
add_test(basis basisTest)

add_executable(systemTest systemTest.cpp)
target_link_libraries(systemTest ${Boost_LIBRARIES})
add_test(system systemTest)

add_executable(eigenHelpersTest eigenHelpersTest.cpp)
target_link_libraries(eigenHelpersTest ${Boost_LIBRARIES})
add_test(eigenHelpers eigenHelpersTest)

add_executable(qcaTest qcaTest.cpp)
target_link_libraries(qcaTest ${Boost_LIBRARIES})
add_test(qca qcaTest)

add_executable(cqcaTest cqcaTest.cpp)
target_link_libraries(cqcaTest ${Boost_LIBRARIES})
add_test(cqca cqcaTest)

add_executable(qcaIsingTest qcaIsingTest.cpp)
target_link_libraries(qcaIsingTest ${Boost_LIBRARIES})
add_test(qcaIsing qcaIsingTest)

add_dependencies(
    buildtest
    utilitiesTest
    basisTest
    systemTest
    eigenHelpersTest
    qcaTest
    cqcaTest
    qcaIsingTest)
