########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(ParallelPerformance ${sources} ${header})

target_link_libraries(ParallelPerformance PRIVATE pz)
if(USING_LOG4CXX)
      target_compile_definitions(ParallelPerformance PRIVATE LOG4CXX)
      target_link_libraries(ParallelPerformance PRIVATE ${Log4cxx_LIBRARY})
 endif()
