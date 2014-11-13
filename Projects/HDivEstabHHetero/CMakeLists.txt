
file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(HDivEstabHHetero ${source} ${headers})

target_link_libraries(HDivEstabHHetero pz)

