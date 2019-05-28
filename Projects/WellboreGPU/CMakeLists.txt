# Minimum cmake version required
cmake_minimum_required (VERSION 3.11.0)

# Project name
project (IntPointsFEM)

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/include /usr/local/include /usr/lib /usr/lib64 /usr/local/lib /usr/local/pzlib/include /usr/local/pzlib/lib .. )

# C++ version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enabling NEOPZ library
IF (WIN32)
    # Prevents timespec redefinition problem with Visual Studio 2015
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
    find_package(PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -mfma -ffast-math -std=c++11") #ligar -O3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
    find_package(PZ REQUIRED)
ENDIF()
include_directories(${PZ_INCLUDE_DIRS})

# Enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread ${PROJECT_SOURCE_DIR}/../externallibs/pthread/include ${PROJECT_SOURCE_DIR}/../pthread )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
    set (PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
    set (PTHREAD_LIB "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})

# Enabling Open MP
option(USING_OMP off)
if(USING_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_OMP")
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif(OPENMP_FOUND)
endif(USING_OMP)

# Enabling CUDA
option(USING_CUDA on)
if(USING_CUDA)
    find_package(CUDA REQUIRED)
    enable_language(CXX CUDA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_CUDA")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --default-stream per-thread  -ftz=true -Xcompiler -fopenmp -O3 -w -Wno-deprecated-gpu-targets")

    set(CMAKE_CUDA_STANDARD 11)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=compute_35 -rdc=true -lcudadevrt)

    set(CUDA_SEPARABLE_COMPILATION ON)
    
    include_directories(/usr/local/cuda-8.0/include/)
    link_directories(/usr/local/cuda-8.0/lib64/)

    file(GLOB sources *.cpp)
    file(GLOB headers *.h)
    file(GLOB cuda *.cu)
    set(main_SRCS ${sources} ${headers} ${cuda})
else()
    file(GLOB sources *.cpp)
    file(GLOB headers *.h)
    set(main_SRCS ${sources} ${headers})
endif(USING_CUDA)

# Configure the path.h.in file.
# This file can only use the @VARIABLE@ syntax in the file
configure_file(path.h.in ${PROJECT_SOURCE_DIR}/path.h @ONLY)

# Create executable and link libraries
if(USING_CUDA)
    cuda_add_executable(IntPointsFEM ${main_SRCS})
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cublas_device_LIBRARY} ${CUDA_cudadevrt_LIBRARY})
elseif(USING_OMP)
    add_executable(IntPointsFEM ${main_SRCS})
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB} OpenMP::OpenMP_CXX)
else()
    add_executable(IntPointsFEM ${main_SRCS})
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB})
endif()

