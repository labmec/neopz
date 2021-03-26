#
# This function will check if any linking is needed when
# using functions from <atomic>
# Usage:
#     check_for_atomic(target)
#
include(CMakePushCheckState)
include(CheckCXXSourceCompiles)

function(check_for_atomic target)
cmake_push_check_state()

if(CMAKE_CXX17_STANDARD_COMPILE_OPTION)
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_CXX17_STANDARD_COMPILE_OPTION})
endif()

check_cxx_source_compiles("
#include <atomic>
std::atomic<int64_t> x;

int main() {
  x.load(std::memory_order_relaxed);
  return std::atomic_is_lock_free(&x);
}"
  atomic)
if(NOT atomic)
  message(STATUS "Failed to compile dummy program with std::atomic")
  message(STATUS "Not a problem: trying to find a lib...")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll" ".so" ".so.1")
  find_library(lib_atomic NAMES atomic)
  if(lib_atomic)
    message(STATUS "Found lib atomic: ${lib_atomic}")
    target_link_libraries(${target} PRIVATE ${lib_atomic})
  else()
    message(FATAL_ERROR "Could not find atomic lib")
  endif()
endif()
cmake_pop_check_state()
mark_as_advanced(lib_atomic)
endfunction()
