
set(main_SRCS QuadraticGeometricTests.cpp)

add_executable(QuadraticGeometricTests ${main_SRCS})

target_link_libraries(QuadraticGeometricTests pz pthread)

