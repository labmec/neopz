## - Config complementary file for the PZ package
# This file only installed (and therefore included by external project) if PZ package has been
# generated with Debug configuration. It WILL NOT BE COPIED for Release configuration of PZ package.
# It adds a definition "PZDEBUG" for CXX Debug configuration. 

if (NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES "-DPZDEBUG")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPZDEBUG" CACHE STRING "Flags used by the CXX compiler during DEBUG builds." FORCE)
endif()