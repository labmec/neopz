include(FindPackageHandleStandardArgs)

find_path(PAPI_ROOT_DIR
  NAMES include/papi.h
  PATHS ENV PAPI_ROOT ${EXTRA_SEARCH_DIRS}
  DOC "PAPI root dir")

find_path(PAPI_INCLUDE_DIR
  NAMES papi.h
  HINTS ${PAPI_ROOT_DIR}
  PATH_SUFFIXES include
  DOC "PAPI include dir")

find_library(PAPI_LIB
  NAMES papi
  HINTS ${PAPI_ROOT_DIR}
  DOC "PAPI library")

mark_as_advanced(PAPI_ROOT_DIR PAPI_INCLUDE_DIR PAPI_LIB)
