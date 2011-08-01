file(GLOB source *.cpp)

set(pzclass5 ${source} PARENT_SCOPE)

add_executable(Class_5 ${source})

target_link_libraries(Class_5 pz pthread)

