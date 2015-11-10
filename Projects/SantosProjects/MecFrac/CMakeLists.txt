#Project test for adaptivity

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(MecFrac ${main_SRCS})

target_link_libraries(MecFrac pz)