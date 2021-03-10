include(FindPackageHandleStandardArgs)

find_path(NUMA_ROOT_DIR
  NAMES include/numa.h
  PATHS ENV NUMA_ROOT ${EXTRA_SEARCH_DIRS}
  DOC "NUMA root dir")

find_path(NUMA_INCLUDE_DIR
  NAMES numa.h
  HINTS ${NUMA_ROOT_DIR}
  PATH_SUFFIXES include
  DOC "NUMA include dir")

find_library(NUMA_LIB
  NAMES numa
  HINTS ${NUMA_ROOT_DIR}
  DOC "NUMA library")

mark_as_advanced(NUMA_ROOT_DIR NUMA_LIB NUMA_INCLUDE_DIR)
