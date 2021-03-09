function(enable_log4cxx target)
    target_link_libraries(${target} PRIVATE ${Log4cxx_LIBRARY})
    target_include_directories(${target} PRIVATE ${Log4cxx_INCLUDE_DIR})
    # if(WIN32) ##TODOWIN32: is it still needed?
    # target_link_libraries(pz PUBLIC odbc32.lib ws2_32.lib mswsock.lib)
    # endif()
endfunction()
