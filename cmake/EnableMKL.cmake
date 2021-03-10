function(enable_mkl target)
  set(CMAKE_MKL_DEBUG OFF)
  find_package(MKL REQUIRED)
  target_link_libraries(${target} PRIVATE mkl::mkl_intel_32bit_omp_dyn)
  target_include_directories(${target} PRIVATE ${MKL_INCLUDE})
  target_compile_definitions(${target} PRIVATE MKLBLAS)
  target_compile_definitions(${target} PRIVATE MKLLAPACK)
  #TODOWIN32: should we do something with mkl_rt on windows?
endfunction()
