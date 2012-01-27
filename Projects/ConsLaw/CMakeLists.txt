
file(GLOB header *.h)

add_executable(ConsLaw main.cpp ${header})

target_link_libraries(ConsLaw pz pthread)
