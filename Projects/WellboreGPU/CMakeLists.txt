option(USING_CUDA off)

if(USING_CUDA)
    find_package(CUDA REQUIRED)
    enable_language(CXX CUDA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_CUDA")
    set(CMAKE_CUDA_FLAGS "--default-stream per-thread -ftz=true -Xcompiler -fopenmp -O3 -w -Wno-deprecated-gpu-targets")

    set(CMAKE_CUDA_STANDARD 11)
    set(CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)	
    file(GLOB sources _SourcesCUDA/*.cu)
    file(GLOB main *.cu)
else()
    	file(GLOB sources *.cpp)
	file(GLOB headers *.h)
endif(USING_CUDA)

set(main_SRCS ${sources} ${headers})

add_executable(WellboreGPU ${main_SRCS})

if(USING_CUDA)
	target_link_libraries(WellboreGPU pz ${PTHREAD_LIB} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cudadevrt_LIBRARY})
else()
	target_link_libraries(WellboreGPU pz ${PTHREAD_LIB})
	configure_file(gmsh/wellbore.msh wellbore.msh COPYONLY)
endif(USING_CUDA)

