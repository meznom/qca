### generic ###
project(ExactDiagonalization_src)
cmake_minimum_required(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

# get program version from git
add_custom_command(OUTPUT version.hpp
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND bash updateVersion.hpp.sh
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

get_target_property(qca_filename _qca LOCATION)
if (APPLE)
    # This makes sure we always link against the absolute filename of the
    # Python and Boost library. Specifically, I needed this to link properly
    # against the Anaconda Python distribution on OS X.
    foreach(lib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
        execute_process(
            COMMAND 
            bash -c "otool -D ${lib} | tail -n1 | tr -d '\n'"
            OUTPUT_VARIABLE lib_install_name)
        add_custom_command(
            TARGET _qca
            POST_BUILD
            COMMAND install_name_tool -change ${lib_install_name} ${lib} ${qca_filename})
    endforeach(lib)
endif(APPLE)

# The Python module always needs to have the name _qca.so.
add_custom_command(TARGET _qca
    POST_BUILD
    COMMAND mv ${qca_filename} ${CMAKE_SOURCE_DIR}/qca/_qca.so)
