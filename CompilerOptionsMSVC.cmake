#Visual Studio 2015 compiler is extremely verbose regarding warnings
#so we cut down warning level to W0
 
# remove default warning level from CMAKE_CXX_FLAGS_INIT
# According to CMake's documentation, MUST BE DONE BEFORE 'project(PZ)' DECLARATION!!
if (MSVC)
	string (REGEX REPLACE "/W3" "/W0" CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT}")
	set(CMAKE_STATIC_LINKER_FLAGS_INIT "${CMAKE_STATIC_LINKER_FLAGS_INIT} /IGNORE:4221,4006")
	set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} /IGNORE:4098")
endif()
