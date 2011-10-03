########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(PerkinsGonzalez ${sources} ${header})

target_link_libraries(PerkinsGonzalez pz pthread)

