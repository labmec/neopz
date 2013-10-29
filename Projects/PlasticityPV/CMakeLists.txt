file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(PlasticityPVMC ${main_SRCS})

target_link_libraries(PlasticityPVMC pz)