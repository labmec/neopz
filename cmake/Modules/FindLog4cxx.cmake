# Taken from https://github.com/Kitware/vibrant/blob/master/CMake/FindLog4cxx.cmake
# Copyright 2010 by Kitware, Inc. All Rights Reserved. Please refer to
# KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
# Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.

# Locate the system installed Log4cxx
# The following variables will be set:
#
# Log4cxx_FOUND       - Set to true if Log4cxx can be found
# Log4cxx_INCLUDE_DIR - The path to the Log4cxx header files
# Log4cxx_LIBRARY     - The full path to the Log4cxx library

if( Log4cxx_DIR )
  find_package( Log4cxx NO_MODULE )
elseif( NOT Log4cxx_FOUND )
  message(STATUS "Searching for log4cxx/logger.h")
  find_path( Log4cxx_INCLUDE_DIR
    NAMES log4cxx/logger.h
    PATHS ${EXTRA_SEARCH_DIRS})

  message(STATUS "Searching for libLog4cxx")
  find_library( Log4cxx_LIBRARY log4cxx )

  include( FindPackageHandleStandardArgs )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( Log4cxx Log4cxx_INCLUDE_DIR Log4cxx_LIBRARY )
  if( LOG4CXX_FOUND )
    set( Log4cxx_FOUND TRUE )

    # From 0.13 onwards, log4cxx has changed the signature of
    # log4cxx::spi::LocationInfo
    # we will test for that
    set(LOG4CXX_TEST_CPP
      "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/log4cxx_quick_test.cpp")
      file(WRITE ${LOG4CXX_TEST_CPP} "
#include \"log4cxx/logger.h\"
#include \"log4cxx/basicconfigurator.h\"
#include \"log4cxx/propertyconfigurator.h\"

int main()
{
  constexpr char fileName[]{\"fileName\"};
  constexpr char functionName[]{\"functionName\"};
  constexpr char shortFunctionName[]{\"shortFunctionName\"};
  constexpr int lineNumber{0};
  log4cxx::spi::LocationInfo(fileName,functionName, shortFunctionName,lineNumber);
  return 0;
}
")
      #so it wont try to link
      set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
      
      # Try to run test program
      try_compile(
        LOG4CXX_TEST_COMPILED
        ${CMAKE_CURRENT_BINARY_DIR}
        ${LOG4CXX_TEST_CPP}
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${Log4cxx_INCLUDE_DIR}"
        OUTPUT_VARIABLE LOG4CXX_TEST_COMPILE_OUTPUT)
    # Check program output
    if (LOG4CXX_TEST_COMPILED)
      set(Log4cxx_NEWER_VERSION TRUE)
    else()
      set(LOG4CXX_TEST_CPP
        "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/log4cxx_quick_test.cpp")
      file(WRITE ${LOG4CXX_TEST_CPP} "
#include \"log4cxx/logger.h\"
#include \"log4cxx/basicconfigurator.h\"
#include \"log4cxx/propertyconfigurator.h\"

int main()
{
  constexpr char fileName[]{\"fileName\"};
  constexpr char functionName[]{\"functionName\"};
  constexpr int lineNumber{0};
  log4cxx::spi::LocationInfo(fileName, functionName,lineNumber);
  return 0;
}
")

      
      # Try to run test program
      try_compile(
        LOG4CXX_TEST_COMPILED
        ${CMAKE_CURRENT_BINARY_DIR}
        ${LOG4CXX_TEST_CPP}
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${Log4cxx_INCLUDE_DIR}"
        OUTPUT_VARIABLE LOG4CXX_TEST_COMPILE_OUTPUT)
      if (LOG4CXX_TEST_COMPILED)
        set(Log4cxx_NEWER_VERSION FALSE)
        message(STATUS "OLDER LOG")
      else()
        message(STATUS "Could not compile test program with log4cxx. Proceeding anyway")
      endif()
    endif()
  endif()
endif()


mark_as_advanced(Log4cxx_INCLUDE_DIR Log4cxx_LIBRARY Log4cxx_NEWER_VERSION)
