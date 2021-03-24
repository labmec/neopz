function(enable_sphinx target)
    
    find_package(Sphinx REQUIRED)
    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/docs_sphinx)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs_sphinx)
    message(${DOXYGEN_OUTPUT_DIR})
    add_custom_target(sphinx
        COMMAND
        ${SPHINX_EXECUTABLE} -E -b html
        -Dbreathe_projects.NeoPZ=${DOXYGEN_OUTPUT_DIR}/xml
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")
endfunction()