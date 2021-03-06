file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(GFlopsTest ${source} ${header})
install(TARGETS GFlopsTest DESTINATION bin)
target_link_libraries(GFlopsTest pz)

