########### next target ###############
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(ValidacaoHCurlFran2 ${main_SRCS})
target_link_libraries(ValidacaoHCurlFran2 pz)