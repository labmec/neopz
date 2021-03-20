function(enable_mkl target)
  #this is for debugging any problems while finding mkl
  set(CMAKE_MKL_DEBUG OFF)
  find_package(MKL REQUIRED)

  set(MKL_THREAD_MODEL_TYPES "Options are: seq, tbb or omp(default)"
    "\nseq: Sequential\ntbb:TBB\nOMP:OpenMP")
  if(NOT MKL_THREAD_MODEL AND TARGET mkl::mkl_intel_32bit_tbb_dyn)
    set(MKL_THREAD_MODEL "tbb" CACHE STRING "${MKL_THREAD_MODEL_TYPES}")
  else()
    set(MKL_THREAD_MODEL "omp" CACHE STRING "${MKL_THREAD_MODEL_TYPES}")
    endif()
  
  if(NOT MKL_THREAD_MODEL STREQUAL "tbb" AND
      NOT MKL_THREAD_MODEL STREQUAL "omp" AND
      NOT MKL_THREAD_MODEL STREQUAL "seq")
    message(FATAL_ERROR "MKL THREADING MODEL ERROR"
      "\n${MKL_THREAD_MODEL_TYPES}")
  endif()
  message(STATUS "Setting MKL threading model: ${MKL_THREAD_MODEL}")
  
  #do NOT change this lib unless you know what you are doing
  target_link_libraries(${target} PRIVATE mkl::mkl_intel_32bit_${MKL_THREAD_MODEL}_dyn)
  target_include_directories(${target} PRIVATE ${MKL_INCLUDE_DIR})
  target_compile_definitions(${target} PRIVATE USING_MKL)
  target_compile_definitions(${target} PRIVATE USING_BLAS)
  target_compile_definitions(${target} PRIVATE USING_LAPACK)
  target_compile_definitions(${target} INTERFACE PZ_USING_MKL)
  target_compile_definitions(${target} INTERFACE PZ_USING_BLAS)
  target_compile_definitions(${target} INTERFACE PZ_USING_LAPACK)
  target_compile_definitions(${target} PRIVATE MKLBLAS)
  target_compile_definitions(${target} PRIVATE MKLLAPACK)
  #TODOWIN32: should we do something with mkl_rt on windows?
  mark_as_advanced(MKL_THREAD_MODEL)
endfunction()