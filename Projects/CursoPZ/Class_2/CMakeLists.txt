file(GLOB source *.cpp)

set(pzclass2 ${source} PARENT_SCOPE)

add_executable(Class_2 ${source})

target_link_libraries(Class_2 pz pthread)

