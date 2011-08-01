file(GLOB source *.cpp)

set(pzclass4 ${source} PARENT_SCOPE)

add_executable(Class_4 ${source})

target_link_libraries(Class_4 pz pthread)

