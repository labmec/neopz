
file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(GeoBlendTests ${source} ${header})

target_link_libraries(GeoBlendTests pz)