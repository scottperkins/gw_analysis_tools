cmake_minimum_required (VERSION 3.11)
project(gw_analysis_tools
	VERSION 0.0 	
	DESCRIPTION "Tools to help with gravitational wave analysis and prediction"
	LANGUAGES CXX C)

# The version number.
set (GWAT_VERSION_MAJOR 1)
set (GWAT_VERSION_MINOR 0)

set (GWAT_ROOT_DIRECTORY "${PROJECT_SOURCE_DIR}/")

option(TEST_PROGRAMS "test-programs" ON)
option(DOCS "DOCS" OFF)
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/gwat/GWATConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/gwat/GWATConfig.h"
  )

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	#Only compile if this project is the main project
	#Just leave empty for now
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

#Required packages
find_package(HDF5 REQUIRED)
find_package(ADOLC REQUIRED)
find_package(OpenMP REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(GSL REQUIRED)

#Library
add_subdirectory("${CMAKE_SOURCE_DIR}/src")

#Executables
add_subdirectory("${CMAKE_SOURCE_DIR}/apps")

if(TEST_PROGRAMS)
	add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
endif()


find_package(Doxygen)
if(DOCS AND DOXYGEN_FOUND)
	add_subdirectory(docs)	
endif()
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
#include_directories("${PROJECT_BINARY_DIR}")
