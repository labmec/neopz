# @file neopz/Projects/RationalFunctionPyramid/CMakeLists.txt  -- Simple project to test implementation of rational functions on pyramidal finite element 

file(GLOB header *.h)

file(GLOB source *.cpp)

#set(pzclass6 ${source} PARENT_SCOPE)

add_executable(RationalFunc ${source} ${header})

target_link_libraries(RationalFunc pz)

