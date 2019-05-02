cmake_minimum_required (VERSION 3.11.0)

project (IntPointsFEM)

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/include /usr/local/include /usr/lib /usr/lib64 /usr/local/lib /usr/local/pzlib/include /usr/local/pzlib/lib .. )

#enabling NEOPZ library
IF (WIN32)
    # Prevents timespec redefinition problem with Visual Studio 2015
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11") #ligar -O3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
ENDIF()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (WIN32)
    find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
    find_package( PZ REQUIRED )
ENDIF (WIN32)


include_directories(${PZ_INCLUDE_DIRS})

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread ${PROJECT_SOURCE_DIR}/../externallibs/pthread/include ${PROJECT_SOURCE_DIR}/../pthread )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
    set (PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
    set (PTHREAD_LIB "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})


#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" ON)
if(USING_MKL)
    ####add_definitions(-DUSING_MKL)
    ####add_definitions(-DUSING_LAPACK)
    #set(USING_LAPACK ON)
    ####add_definitions(-DUSING_BLAS)
    #set(USING_BLAS ON)
    if (WIN32)
        set(MKL_ROOT “opt/intel/compilers_and_libraries/mac”)
        find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
        find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
    else()
        find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
        find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
                ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
        find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
                /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
        find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
                ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
        find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
                /opt/intel/composer_xe/compiler/lib
                /opt/intel/composerxe/lib/intel64
                /opt/intel/lib /softwares/intel/lib/intel64
                /opt/intel/compilers_and_libraries/linux/lib/intel64_lin )

        #find_path(OMP_INCLUDE NAMES omp.h PATHS ${SEARCH_DIRS} /opt/intel/composer_xe_2013_sp1.1.103/compiler/include)
        #include_directories(${OMP_INCLUDE})
    endif()

    if(MKL_INCLUDE-NOTFOUND)
        set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
    else()
        #include_directories(${MKL_INCLUDE})
    endif()
    if(MKL_LIB-NOTFOUND)
        set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
    else()
        link_libraries(${MKL_LIB_INTEL})
        link_libraries(${MKL_LIB_CORE})
        link_libraries(${MKL_LIB_THREAD})
        link_libraries(${COMPOSER_OMP})
    endif()
endif(USING_MKL)

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_TBB")
    if(WIN32)
        set(TBB_ROOT "C:/Program Files (x86)/IntelSWTools/TBB")
        find_path(TBB_INCLUDE tbb.h HINTS ${TBB_ROOT}/tbb2018_20180618oss/include/tbb NO_DEFAULT_PATH)
        find_library(TBB_LIB tbb.lib HINTS ${TBB_ROOT}/tbb2018_20180618oss/lib/intel64/vc14 NO_DEFAULT_PATH)

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
        find_file(TBB_HEADER_FILE tbb.h HINTS /usr/local/TBB/tbb2018_20180618oss/include/tbb/ NO_DEFAULT_PATH)
        #--path do include do tbb no cluster
        #find_file(TBB_HEADER_FILE tbb.h /softwares/intel/tbb/include/tbb NO_DEFAULT_PATH)
        #-----------------------------------
        if(NOT TBB_HEADER_FILE MATCHES "tbb.h")
            set (TBB_HEADER_FILE "" CACHE PATH "Path to 'tbb.h' file")
            message(FATAL_ERROR "TBB_HEADER_FILE must contain a valid path to 'tbb.h' file!")
        else()
            get_filename_component(TBB_ROOT_DIR ${TBB_HEADER_FILE} DIRECTORY)
            get_filename_component(TBB_INCLUDE ${TBB_ROOT_DIR} DIRECTORY)
            get_filename_component(TBB_ROOT_DIR ${TBB_INCLUDE} DIRECTORY)
            if(APPLE)
                set(TBB_LIB ${TBB_ROOT_DIR}/lib/libtbb.dylib)
            else() # Linux
                set(TBB_LIB ${TBB_ROOT_DIR}/lib/intel64/gcc4.7/libtbb.so)
                #--path da lib do tbb no cluster
                #set(TBB_LIB ${TBB_ROOT_DIR}/lib/intel64_lin/gcc4.4/libtbb.so)
                #-------------------------------
            endif()

            include_directories(${TBB_INCLUDE})
            link_libraries(${TBB_LIB})
        endif()
    endif(WIN32)

endif(USING_TBB)

#Specify where the headers files must be found
include_directories("${PROJECT_SOURCE_DIR}/_Headers")
set(PROJECTSOURCE ${BASIC_TYPES} ${HEADERS})


file(GLOB headers *.h)
file(GLOB sources *.cpp)

#enabling cuda
option(USING_CUDA on)
if(USING_CUDA)
    find_package(CUDA REQUIRED)
    enable_language(CXX CUDA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_CUDA")
    set(CMAKE_CUDA_FLAGS "--default-stream per-thread -ftz=true -Xcompiler -fopenmp -O3 -w -Wno-deprecated-gpu-targets")

    set(CMAKE_CUDA_STANDARD 11)
    set(CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)	

    file(GLOB cuda *.cu)
    list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/TPZIntPointsFEM.cpp)
    set(main_SRCS ${sources} ${headers} ${cuda})
else()
    set(main_SRCS ${sources} ${headers})
endif(USING_CUDA)

# Configure the path.h.in file.
# This file can only use the @VARIABLE@ syntax in the file
configure_file(path.h.in ${PROJECT_SOURCE_DIR}/path.h @ONLY)

add_executable(IntPointsFEM ${main_SRCS})

if(USING_CUDA)
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY})
    configure_file(gmsh/wellbore.msh wellbore.msh COPYONLY)
    configure_file(gmsh/wellbore-3.msh wellbore-3.msh COPYONLY)
else()
    target_link_libraries(IntPointsFEM pz ${PTHREAD_LIB})
endif(USING_CUDA)

