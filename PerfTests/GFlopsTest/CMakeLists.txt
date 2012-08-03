file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(GFlopsTest ${source} ${header})

target_link_libraries(GFlopsTest pz pthread)

