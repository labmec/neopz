########### next target ###############

file(GLOB sources *.cpp)

set(main_SRCS ${sources})

add_executable(PerkinsGonzalez ${main_SRCS})

target_link_libraries(PerkinsGonzalez pz pthread)

