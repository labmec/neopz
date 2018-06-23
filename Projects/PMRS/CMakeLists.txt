
# Adding projec files

file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(PMRS ${source} ${header})

target_link_libraries(PMRS pz)




