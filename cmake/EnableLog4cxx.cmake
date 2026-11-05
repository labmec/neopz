function(enable_log4cxx target)
  find_package(Log4cxx REQUIRED)
  target_link_libraries(${target} PRIVATE ${Log4cxx_LIBRARY})
  target_include_directories(${target} PRIVATE ${Log4cxx_INCLUDE_DIR})
  target_compile_definitions(${target} PUBLIC PZ_LOG)
  if(Log4cxx_NEWER_VERSION)
    target_compile_definitions(${target} PRIVATE PZ_LOG_LOG4CXXNEW)
    message(STATUS "Found Log4cxx newer than 0.13 (or equal)")
  else()
    target_compile_definitions(${target} PRIVATE PZ_LOG_LOG4CXXOLD)
  endif()
  
  
  # if(WIN32) ##TODOWIN32: is it still needed?
  # target_link_libraries(pz PUBLIC odbc32.lib ws2_32.lib mswsock.lib)
  # endif()
endfunction()
