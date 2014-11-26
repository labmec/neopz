########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(GlobalPostPro ${sources} ${header})

target_link_libraries(GlobalPostPro pz)
