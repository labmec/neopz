
file(GLOB header *.h)

add_executable(EulerTest main.cpp ${header})

target_link_libraries(EulerTest pz pthread)
