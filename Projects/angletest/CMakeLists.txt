
set(testeangular_SRCS testeangular.cpp)

add_executable(testeangular ${testeangular_SRCS})

target_link_libraries(testeangular pz pthread)

