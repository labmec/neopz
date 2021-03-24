file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(MatLib ${Mymain})

target_link_libraries(MatLib pz ${MATLIB_LIB})

