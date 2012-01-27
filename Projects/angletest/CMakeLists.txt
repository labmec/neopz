
file(GLOB header *.h)

set(testeangular_SRCS testeangular.cpp)

add_executable(angletest ${testeangular_SRCS} ${header})

target_link_libraries(angletest pz pthread)

