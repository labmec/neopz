file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(adaptivity ${sources} ${header})

target_link_libraries(adaptivity pz pthread)

