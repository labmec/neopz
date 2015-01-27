########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Itasca ${sources} ${header})

target_link_libraries(Itasca pz) 

file(COPY example1.dat DESTINATION .)
file(COPY "example1 BIG.dat" DESTINATION .)

