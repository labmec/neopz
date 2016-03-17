set(main_SRCS testMatLab.C)

add_executable(testMatLab ${main_SRCS})

target_link_libraries(testMatLab pz)