file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(ExemploMapGeo ${main_SRCS})

target_link_libraries(ExemploMapGeo pz)