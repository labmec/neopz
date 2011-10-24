set(main_SRCS Cubo_do_Nathan.cpp pzpostprocanalysis.h pzpostprocanalysis.cpp pzcompelpostproc.h pzpostprocmat.h pzpostprocmat.cpp)

add_executable(Cubo_do_Nathan ${main_SRCS})

target_link_libraries(Cubo_do_Nathan pz pthread)

file(COPY cube1.txt DESTINATION .)

