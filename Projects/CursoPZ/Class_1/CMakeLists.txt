file(GLOB source *.cpp)

set(pzclass1 ${source} PARENT_SCOPE)

add_executable(Class_1 ${source})

target_link_libraries(Class_1 pz pthread)

