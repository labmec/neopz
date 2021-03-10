include(FindPackageHandleStandardArgs)

find_path(LIKWID_ROOT_DIR
  NAMES include/likwid.h
  PATHS ENV LIKWID_ROOT ${EXTRA_SEARCH_DIRS}
  DOC "LIKWID root dir")

find_path(LIKWID_INCLUDE_DIR
  NAMES likwid.h
  HINTS ${LIKWID_ROOT_DIR}
  PATH_SUFFIXES include
  DOC "LIKWID include dir")

find_library(LIKWID_LIB
  NAMES likwid
  HINTS ${LIKWID_ROOT_DIR}
  DOC "LIKWID library")

mark_as_advanced(LIKWID_ROOT_DIR LIKWID_INCLUDE_DIR LIKWID_LIB)
