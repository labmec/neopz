# @file neopz/Projects/descontinuo/CMakeLists.txt  -- CMake file for the discontinuous approximation function tutorial 

file(GLOB source *.cpp)
file(GLOB header *.h)

set(discontinuous ${discontinuous} ${source} ${header} PARENT_SCOPE)

# First executable testing Euler equation
add_executable(discontinuous_Euler euler.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

target_link_libraries(discontinuous_Euler pz pthread)

# Second executable testing Euler equation
#add_executable(discontinuous_Euler0 euler0.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

#target_link_libraries(discontinuous_Euler0 pz pthread)

#Executable to three dimensional test
add_executable(discontinuous_Wing wingnaca0012.cpp TPZFlowCMesh.cpp TPZIterativeAnalysis.cpp ${header})

target_link_libraries(discontinuous_Wing pz pthread)
