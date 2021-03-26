file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain MatrixMarket.cpp ${sources} ${headers})

add_executable(MatrixMarket ${Mymain})

target_link_libraries(MatrixMarket pz)

