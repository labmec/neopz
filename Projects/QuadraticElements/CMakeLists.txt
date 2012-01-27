
file(GLOB header *.h)

set(main_SRCS QuadraticGeometricTests.cpp)

add_executable(QuadraticElements ${main_SRCS} ${header})

target_link_libraries(QuadraticElements pz pthread)

