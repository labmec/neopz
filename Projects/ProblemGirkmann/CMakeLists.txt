file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(ProblemGirkmann ${main_SRCS} ${headers})

target_link_libraries(ProblemGirkmann pz pthread)