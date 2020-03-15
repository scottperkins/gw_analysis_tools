# - Try to find Adolc
# Once done this will define
#  FFTW3_FOUND - System has Adolc
#  FFTW3_INCLUDE_DIRS - The Adolc include directories
#  FFTW3_LIBRARY_DIRS - The library directories needed to use Adolc
#  FFTW3_LIBRARIES    - The libraries needed to use Adolc
#  Adapted from https://raw.githubusercontent.com/joaoleal/CppADCodeGen/master/cmake/FindADOLC.cmake

IF (FFTW3_INCLUDES AND FFTW3_LIBRARIES)
  SET(FFTW3_FIND_QUIETLY TRUE FFTW3)
ENDIF (FFTW3_INCLUDES AND FFTW3_LIBRARIES)


FIND_PACKAGE(PkgConfig)
IF( PKG_CONFIG_FOUND )
  pkg_check_modules( FFTW3 QUIET fftw3>=3)
ENDIF( PKG_CONFIG_FOUND )


IF( FFTW3_FOUND )
  IF(NOT FFTW3_FIND_QUIETLY)
    MESSAGE(STATUS "package FFTW3 ${FFTW3_VERSION} found")
  ENDIF()
ELSE( FFTW3_FOUND )
  string(REPLACE ":" ";" INC_SEARCH_LIST $ENV{CPATH})
  FIND_PATH(FFTW3_INCLUDE_DIRS NAMES fftw3.h
            HINTS  
                   "/usr/include" 
                   "/usr/local/include" 
		   ${INC_SEARCH_LIST}
		   )
           
  string(REPLACE ":" ";" LIB_SEARCH_LIST $ENV{LIBRARY_PATH})
  FIND_LIBRARY(FFTW3_LIBRARY fftw3
                HINTS 
                      "/usr/lib" 
                      "/usr/local/lib" 
		      ${LIB_SEARCH_LIST})
  #CHECK_SYMBOL_EXISTS(function ADOLC_LIBRARY FUNCTION_FOUND)
   
  SET(FFTW3_INCLUDE_DIR ${FFTW3_INCLUDE_DIRS})
  SET(FFTW3_LIBRARIES ${FFTW3_LIBRARY})

  INCLUDE(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(FFTW3  DEFAULT_MSG
                                    FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)

  MARK_AS_ADVANCED(FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)
  
  IF( FFTW3_FOUND AND NOT FFTW3_FIND_QUIETLY )
    MESSAGE(STATUS "package fftw3 found")
  ENDIF()
ENDIF( FFTW3_FOUND )
