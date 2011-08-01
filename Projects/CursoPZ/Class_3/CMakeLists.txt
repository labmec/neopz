file(GLOB source *.cpp)

set(pzclass3 ${source} PARENT_SCOPE)

add_executable(Class_3 ${source})

target_link_libraries(Class_3 pz pthread)

