function(enable_sphinx target)
    find_package(Sphinx REQUIRED)
    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

    set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in)
    set(SPHINX_OUT ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

    set(REFS_IN ${CMAKE_CURRENT_SOURCE_DIR}/references.bib)
    set(REFS_OUT ${CMAKE_CURRENT_BINARY_DIR}/references.bib)
    configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)
    configure_file(${REFS_IN} ${REFS_OUT} @ONLY)
    add_custom_target(run-sphinx
        COMMAND
        ${SPHINX_EXECUTABLE} -E -b html -c ${CMAKE_CURRENT_BINARY_DIR}
        -Dbreathe_projects.NeoPZ=${DOXYGEN_OUTPUT_DIR}/xml
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")
endfunction()