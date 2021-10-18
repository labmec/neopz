function(enable_doxygen target)
    if (NOT "$ENV{DOXYGEN_ROOT}" STREQUAL "")
        message(STATUS "Looking for doxygen in $ENV{BOOST_ROOT}")
    endif()
    find_package(Doxygen
        REQUIRED)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
    if(USING_LOG4CXX)
        set(DOX_PZ_LOG 1)
    endif()
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(run-doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )

    ##so sphinx can find it
    set(DOXYGEN_OUTPUT_DIR ${DOXYGEN_OUTPUT_DIR} PARENT_SCOPE)
endfunction()
