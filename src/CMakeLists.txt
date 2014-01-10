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
if (APPLE)
    # This makes sure we always link against the absolute filename of the
    # Python library. Speficially, I needed this to link properly against the
    # Anaconda Python distribution on OS X.
    execute_process(COMMAND bash -c "otool -D ${PYTHON_LIBRARIES} | tail -n1 | tr -d '\n'" OUTPUT_VARIABLE INSTALL_NAME)
    add_custom_command(TARGET _qca
        POST_BUILD
        COMMAND install_name_tool -change ${INSTALL_NAME} ${PYTHON_LIBRARIES} ${LIB_FILENAME}
        COMMAND mv ${LIB_FILENAME} ${CMAKE_SOURCE_DIR}/qca/_qca.so)
else()
    add_custom_command(TARGET _qca
        POST_BUILD
        COMMAND mv ${LIB_FILENAME} ${CMAKE_SOURCE_DIR}/qca/_qca.so)
endif()
