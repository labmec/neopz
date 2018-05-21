#Visual Studio 2015 compiler is extremely verbose regarding warnings
#so we cut down warning level to W0
 
# remove default warning level from CMAKE_CXX_FLAGS_INIT
# According to CMake's documentation, MUST BE DONE BEFORE 'project(PZ)' DECLARATION!!
if (MSVC)
	string (REGEX REPLACE "/W3" "/W0" CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT}")
endif()
