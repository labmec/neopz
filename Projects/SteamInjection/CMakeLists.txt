file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(SteamInjection ${main_SRCS})

target_link_libraries(SteamInjection pz pthread)

