
file(GLOB header *.h)

set(main_SRCS QuadraticGeometricTests.cpp)

add_executable(QuadraticGeometricTests ${main_SRCS} ${header})

target_link_libraries(QuadraticGeometricTests pz pthread)

