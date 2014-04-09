# @file neopz/Projects/RationalFunction/CMakeLists.txt  -- Simple project to test implementation of rational functions on pyramidal finite element with performance evaluation

file(GLOB header *.h)

file(GLOB source *.cpp)

#set(pzclass6 ${source} PARENT_SCOPE)

add_executable(RationalFunc-Perf ${source} ${header})

target_link_libraries(RationalFunc-Perf pz)

