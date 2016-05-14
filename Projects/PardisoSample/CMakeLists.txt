########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp pardiso_mkl.cpp)

add_executable(PardisoSample ${main_SRCS} ${header})

target_link_libraries(PardisoSample pz ${LIB_PARDISO})