
file(GLOB header *.h)

set(main_SRCS Bima1.cpp)

add_executable(HibridoBima ${main_SRCS} ${header})

target_link_libraries(HibridoBima pz)

