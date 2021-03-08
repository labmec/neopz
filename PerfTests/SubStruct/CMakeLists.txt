file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-SubStruct ${source} ${header})
install(TARGETS Perf-SubStruct DESTINATION bin)
target_link_libraries(Perf-SubStruct PRIVATE pz)

configure_file(8andares02.txt 8andares02.txt COPYONLY)
configure_file(cube1.txt cube1.txt COPYONLY)

if(USING_LOG4CXX)
      target_compile_definitions(Perf-SubStruct PRIVATE LOG4CXX)
      target_link_libraries(Perf-SubStruct PRIVATE ${Log4cxx_LIBRARY})
endif()
