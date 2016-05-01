########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp)

add_executable(PardisoSample ${main_SRCS} ${header})

target_link_libraries(PardisoSample pz ${LIB_PARDISO})