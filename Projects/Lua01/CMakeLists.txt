########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Lua01 ${sources} ${header})

target_link_libraries(Lua01 pz pthread)

