project(diagonaliser_tests)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(..)

add_executable(utilitiesTest utilitiesTest.cpp)
target_link_libraries(utilitiesTest boost_unit_test_framework)
add_test(utilities utilitiesTest)
add_dependencies(buildtest utilitiesTest)
