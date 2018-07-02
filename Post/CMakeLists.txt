# @file neopz/Post/CMakeLists.txt  -- CMake file for the Post module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)

if(BUILD_COMPLEX_PROJECTS)
    file (GLOB to_remove_from_headers pzgradientreconstruction.h)
    list (REMOVE_ITEM headers ${to_remove_from_headers})
endif() 

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Post)

file(GLOB sources *.cpp)

if(BUILD_COMPLEX_PROJECTS)
    file (GLOB to_remove_from_sources pzgradientreconstruction.cpp)
    list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

set(pzpost_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


