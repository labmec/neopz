# @file neopz/Projects/descontinuo/CMakeLists.txt  -- CMake file for the discontinuous approximation function tutorial 

file(GLOB source *.cpp)
file(GLOB header *.h)

set(discontinuous ${discontinuous} ${source} ${header} PARENT_SCOPE)

# First executable testing Euler equation
add_executable(DiscontinuousEuler euler.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

target_link_libraries(DiscontinuousEuler pz pthread)

# Second executable testing Euler equation
#add_executable(DiscontinuousEuler0 euler0.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

#target_link_libraries(DiscontinuousEuler0 pz pthread)

#Executable to three dimensional test
add_executable(DiscontinuousWing wingnaca0012.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

target_link_libraries(DiscontinuousWing pz pthread)
