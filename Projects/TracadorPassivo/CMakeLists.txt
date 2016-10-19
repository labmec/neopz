
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(TracadorPassivo ${source} ${header})

target_link_libraries(TracadorPassivo pz)

