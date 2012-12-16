file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(PoroElasticAnalysis ${Mymain})

target_link_libraries(PoroElasticAnalysis pz)

file(COPY Input DESTINATION .)

