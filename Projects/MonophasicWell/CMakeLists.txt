
file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(MonoPhasicWell ${source} ${header})

target_link_libraries(MonoPhasicWell pz)