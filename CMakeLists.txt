cmake_minimum_required (VERSION 2.6)
project(gw_analysis_tools)

# The version number.
set (GWAT_VERSION_MAJOR 1)
set (GWAT_VERSION_MINOR 0)

set (GWAT_ROOT_DIRECTORY "${PROJECT_SOURCE_DIR}/")
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/docs/GWATConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/GWATConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
#include_directories("${PROJECT_BINARY_DIR}")
