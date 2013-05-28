### generic ###
project(ExactDiagonalization_src)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

# get program version from git
add_custom_command(OUTPUT version.hpp
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND bash updateVersion.hpp.sh
)

# runQca
add_executable(runQca
    runQca.cpp
    version.hpp
)

# Python binding
add_library(_qca
    SHARED
    pythonbinding.cpp
    qca.hpp
    version.hpp)

target_link_libraries(_qca
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES})

# The Python module always needs to have the name _qca.so.
get_target_property(LIB_FILENAME _qca LOCATION)
add_custom_command(TARGET _qca
    POST_BUILD         
    COMMAND mv ${LIB_FILENAME} ${CMAKE_SOURCE_DIR}/qca/_qca.so)
