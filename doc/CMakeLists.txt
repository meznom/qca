project(ExactDiagonalisationDoc)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

add_custom_target(doc
    COMMAND doxygen
)
