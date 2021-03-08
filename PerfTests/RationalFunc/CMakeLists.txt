# @file neopz/Projects/RationalFunction/CMakeLists.txt  -- Simple project to test implementation of rational functions on pyramidal finite element with performance evaluation

file(GLOB header *.h)

file(GLOB source *.cpp)

#set(pzclass6 ${source} PARENT_SCOPE)

add_executable(RationalFunc-Perf ${source} ${header})

target_link_libraries(RationalFunc-Perf PRIVATE pz)
if(USING_LOG4CXX)
      target_compile_definitions(RationalFunc-Perf PRIVATE LOG4CXX)
      target_link_libraries(RationalFunc-Perf PRIVATE ${Log4cxx_LIBRARY})
endif()
