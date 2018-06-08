file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main dfn_poroelastic.cpp ${sources} ${headers})

add_executable(dfn_poroelastic ${main})

target_link_libraries(dfn_poroelastic pz)

