
# Adding projec files

file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(PoroPermCoupling ${source} ${header})

target_link_libraries(PoroPermCoupling pz)




