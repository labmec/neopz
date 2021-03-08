file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-DecomposeLDLt ${source} ${header})
install(TARGETS Perf-DecomposeLDLt DESTINATION bin)
target_link_libraries(Perf-DecomposeLDLt PRIVATE pz)
if(USING_LOG4CXX)
      target_compile_definitions(Perf-DecomposeLDLt PRIVATE LOG4CXX)
      target_link_libraries(Perf-DecomposeLDLt PRIVATE ${Log4cxx_LIBRARY})
endif()
