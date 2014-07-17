file(GLOB source *.cpp)
file(GLOB header *.h)

set(visualmatrix ${visualmatrix} ${source} ${header} PARENT_SCOPE)

add_executable(proDG ${source} ${header})

target_link_libraries(proDG pz)

