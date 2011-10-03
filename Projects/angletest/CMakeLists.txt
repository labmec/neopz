
file(GLOB header *.h)

set(testeangular_SRCS testeangular.cpp)

add_executable(testeangular ${testeangular_SRCS} ${header})

target_link_libraries(testeangular pz pthread)

