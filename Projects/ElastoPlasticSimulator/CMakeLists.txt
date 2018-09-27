
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main main.cpp ${sources} ${headers})

add_executable(ElastoPlasticSimulator ${main})

target_link_libraries(ElastoPlasticSimulator pz)

