
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(PauloBosing ${source} ${header})

target_link_libraries(PauloBosing pz)

