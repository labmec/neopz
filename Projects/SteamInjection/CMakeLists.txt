# @file neopz/Projects/SteamInjection/CMakeLists.txt  -- CMake file for the SteamInjection project 

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(SteamInjection ${main_SRCS} ${headers})

target_link_libraries(SteamInjection pz pthread)

