# @file neopz/External/CMakeLists.txt  -- CMake file for the External module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzexternal_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

add_subdirectory(sloan)
set(sloan_SRCS ${sloan_SRCS} PARENT_SCOPE)

#enabling FAD library
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	#Adding a variable to hold this definition
	add_subdirectory(FAD)
	set(FAD_SRCS ${FAD_SRCS} PARENT_SCOPE)
	set(FAD_fad_SRCS ${FAD_fad_SRCS} PARENT_SCOPE)
	set(FAD_Tiny_SRCS ${FAD_Tiny_SRCS} PARENT_SCOPE)
	set(FAD_TinyET_SRCS ${FAD_TinyET_SRCS} PARENT_SCOPE)
endif(USING_FAD)




