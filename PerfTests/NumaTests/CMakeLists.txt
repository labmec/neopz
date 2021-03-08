file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-NumaTst ${source} ${header})
install(TARGETS Perf-NumaTst DESTINATION bin)
target_link_libraries(Perf-NumaTst PRIVATE pz Threads::Threads)
if(USING_LOG4CXX)
      target_compile_definitions(Perf-NumaTst PRIVATE LOG4CXX)
      target_link_libraries(Perf-NumaTst PRIVATE ${Log4cxx_LIBRARY})
endif()
