project(ExactDiagonalization_tests)
cmake_minimum_required(VERSION 2.6)

include_directories(
	../src
)

find_package(Boost 1.49.0 COMPONENTS unit_test_framework)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -USTORAGE_TYPE_OF_FERMIONIC_STATE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -USTORAGE_TYPE_OF_FERMIONIC_STATE")

set(tests basisTest systemTest qcaTest qcaIsingTest)
foreach(test ${tests})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} ${Boost_LIBRARIES})
    add_test(${test}_test ${test})

    if (APPLE)
        # Use install_name_tool to make sure we link against the absolute path
        # of linked libraries (and not relative paths)
        foreach(lib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
            get_target_property(test_filename ${test} LOCATION)
            execute_process(
                COMMAND
                bash -c "otool -D ${lib} | tail -n1 | tr -d '\n'"
                OUTPUT_VARIABLE lib_install_name)
            add_custom_command(
                TARGET ${test}
                POST_BUILD
                COMMAND install_name_tool -change ${lib_install_name} ${lib} ${test_filename})
        endforeach(lib)
    endif(APPLE)
endforeach(test)

add_dependencies(buildtest ${tests})
