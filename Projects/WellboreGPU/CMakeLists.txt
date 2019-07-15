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
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mfma -ffast-math -std=c++11") #ligar -O3
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

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_TBB")
    if(WIN32)
        set(TBB_ROOT "C:/Program Files (x86)/IntelSWTools/TBB")
        find_path(TBB_INCLUDE tbb.h HINTS ${TBB_ROOT}/tbb2018_20180618oss/include/tbb NO_DEFAULT_PATH)
        find_library(TBB_LIB_RELEASE tbb.lib HINTS ${TBB_ROOT}/tbb2018_20180618oss/lib/intel64/vc14 NO_DEFAULT_PATH)
        find_library(TBB_LIB_DEBUG tbb_debug.lib HINTS ${TBB_ROOT}/tbb2018_20180618oss/lib/intel64/vc14 NO_DEFAULT_PATH)

        set(TBB_LIB
                debug ${TBB_LIB_DEBUG}
                optimized ${TBB_LIB_RELEASE}
                )

        if(TBB_INCLUDE-NOTFOUND)
            set (TBB_INCLUDE "" CACHE PATH "Directory where tbb header files can be found")
        else()
            include_directories(${TBB_INCLUDE}/..)
        endif()
        if(TBB_LIB-NOTFOUND)
            set (TBB_LIB "" CACHE PATH "Directory where the tbb library can be found")
        else()
            link_libraries(${TBB_LIB})
        endif()
    else()
        find_file(TBB_HEADER_FILE tbb.h HINTS /usr/local/TBB/*/include/tbb/ /softwares/intel/tbb/include/tbb/ NO_DEFAULT_PATH)

        if(NOT TBB_HEADER_FILE MATCHES "tbb.h")
            set (TBB_HEADER_FILE "" CACHE PATH "Path to 'tbb.h' file")
            message( FATAL_ERROR "TBB_HEADER_FILE must contain a valid path to 'tbb.h' file!")
        else()
            get_filename_component(TBB_ROOT_DIR ${TBB_HEADER_FILE} DIRECTORY)
            get_filename_component(TBB_INCLUDE ${TBB_ROOT_DIR} DIRECTORY)
            get_filename_component(TBB_ROOT_DIR ${TBB_INCLUDE} DIRECTORY)
            if(APPLE)
                set(TBB_LIB ${TBB_ROOT_DIR}/lib/libtbb.dylib)
            else() # Linux
                find_file(TBB_LIB_FILE libtbb.so HINTS ${TBB_ROOT_DIR}/lib/intel64/gcc4.7/ ${TBB_ROOT_DIR}/lib/intel64_lin/gcc4.4/ NO_DEFAULT_PATH)
                get_filename_component(TBB_LIB_DIR ${TBB_LIB_FILE} DIRECTORY)
                set(TBB_LIB ${TBB_LIB_DIR}/libtbb.so)

            endif()

            include_directories(${TBB_INCLUDE})
            link_libraries(${TBB_LIB})
        endif()
    endif(WIN32)
endif(USING_TBB)

# Residual option
option(USING_SPARSE off)
if(USING_SPARSE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_SPARSE")
endif(USING_SPARSE)

# Assembly options
option(COMPUTE_K_GS off)
if(COMPUTE_K_GS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOMPUTE_K_GS")
endif(COMPUTE_K_GS)

# Polinomial order
option(O_LINEAR on)
if(O_LINEAR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DO_LINEAR")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DO_LINEAR") 
    set(O_QUADRATIC off)
    set(O_CUBIC off)
endif(O_LINEAR)

option(O_QUADRATIC off)
if(O_QUADRATIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DO_QUADRATIC")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DO_QUADRATIC")
    set(O_LINEAR off)
    set(O_CUBIC off)
endif(O_QUADRATIC)

option(O_CUBIC off)
if(O_CUBIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DO_CUBIC")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DO_CUBIC")
    set(O_LINEAR off)
    set(O_QUADRATIC off)
endif(O_CUBIC)

# Enabling CUDA
option(USING_CUDA off)
if(USING_CUDA)
    find_package(CUDA REQUIRED)
    enable_language(CXX CUDA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_CUDA")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --default-stream per-thread  -ftz=true -Xcompiler -O3 --ptxas-options=-O3 -w -Wno-deprecated-gpu-targets -arch=compute_35 -rdc=true -lcudadevrt --use_fast_math")

    set(CMAKE_CUDA_STANDARD 11)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --default-stream per-thread  -ftz=true -Xcompiler -O3 --ptxas-options=-O3 -w -Wno-deprecated-gpu-targets -arch=compute_35 -rdc=true -lcudadevrt --use_fast_math")

    set(CUDA_SEPARABLE_COMPILATION ON)
    
    include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include/)
    link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64/)

    file(GLOB sources *.cpp)
    file(GLOB headers *.h)
    file(GLOB cuda *.cu)
    set(main_SRCS ${sources} ${headers} ${cuda})
else()
    file(GLOB sources *.cpp)
    file(GLOB headers *.h)
    set(main_SRCS ${sources} ${headers})
endif(USING_CUDA)

# K and R assembly using PZ
option(COMPUTE_WITH_PZ off)
if(COMPUTE_WITH_PZ)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOMPUTE_WITH_PZ")
    set(USING_SPARSE off)
    set(COMPUTE_K_GS off)
    set(USING_CUDA off)
endif(COMPUTE_WITH_PZ)

# Configure the path.h.in file.
# This file can only use the @VARIABLE@ syntax in the file
configure_file(path.h.in ${PROJECT_SOURCE_DIR}/path.h @ONLY)

# Create executable and link libraries
if(USING_CUDA)
#    cuda_add_executable(IntPointsFEM ${main_SRCS})
    add_executable(IntPointsFEM ${main_SRCS})
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cublas_device_LIBRARY} ${CUDA_cudadevrt_LIBRARY})
else()
    add_executable(IntPointsFEM ${main_SRCS})
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB})
endif()

