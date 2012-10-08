
set(RefPat_SRCS main.cpp)

add_executable(Test_RefPattern ${RefPat_SRCS})

target_link_libraries(Test_RefPattern pz)

#install(TARGETS Test_RefPattern ${INSTALL_TARGETS_DEFAULT_ARGS})

