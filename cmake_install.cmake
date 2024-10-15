# Install script for directory: /opt/gw_analysis_tools

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/opt/gw_analysis_tools/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gwat" TYPE FILE FILES
    "/opt/gw_analysis_tools/include/gwat/D_Z_Config.h"
    "/opt/gw_analysis_tools/include/gwat/D_Z_Config_modified_dispersion.h"
    "/opt/gw_analysis_tools/include/gwat/EA_IMRPhenomD_NRT.h"
    "/opt/gw_analysis_tools/include/gwat/GWATConfig.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomD.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomD_NRT.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomP.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomP_NRT.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomPv3.h"
    "/opt/gw_analysis_tools/include/gwat/IMRPhenomPv3utils.h"
    "/opt/gw_analysis_tools/include/gwat/QNM_data.h"
    "/opt/gw_analysis_tools/include/gwat/TaylorT2.h"
    "/opt/gw_analysis_tools/include/gwat/autocorrelation.h"
    "/opt/gw_analysis_tools/include/gwat/autocorrelation_cuda.h"
    "/opt/gw_analysis_tools/include/gwat/detector_util.h"
    "/opt/gw_analysis_tools/include/gwat/error.h"
    "/opt/gw_analysis_tools/include/gwat/fisher.h"
    "/opt/gw_analysis_tools/include/gwat/gIMRPhenomD.h"
    "/opt/gw_analysis_tools/include/gwat/gIMRPhenomP.h"
    "/opt/gw_analysis_tools/include/gwat/gwatpy_wrapping.h"
    "/opt/gw_analysis_tools/include/gwat/io_util.h"
    "/opt/gw_analysis_tools/include/gwat/mc_reject.h"
    "/opt/gw_analysis_tools/include/gwat/mcmc_gw.h"
    "/opt/gw_analysis_tools/include/gwat/mcmc_gw_extended.h"
    "/opt/gw_analysis_tools/include/gwat/mcmc_io_util.h"
    "/opt/gw_analysis_tools/include/gwat/mcmc_sampler.h"
    "/opt/gw_analysis_tools/include/gwat/mcmc_sampler_internals.h"
    "/opt/gw_analysis_tools/include/gwat/ortho_basis.h"
    "/opt/gw_analysis_tools/include/gwat/pn_waveform_util.h"
    "/opt/gw_analysis_tools/include/gwat/ppE_IMRPhenomD.h"
    "/opt/gw_analysis_tools/include/gwat/ppE_IMRPhenomD_NRT.h"
    "/opt/gw_analysis_tools/include/gwat/ppE_IMRPhenomP.h"
    "/opt/gw_analysis_tools/include/gwat/ppE_utilities.h"
    "/opt/gw_analysis_tools/include/gwat/standardPriorLibrary.h"
    "/opt/gw_analysis_tools/include/gwat/threadPool.h"
    "/opt/gw_analysis_tools/include/gwat/util.h"
    "/opt/gw_analysis_tools/include/gwat/waveform_generator.h"
    "/opt/gw_analysis_tools/include/gwat/waveform_generator_C.h"
    "/opt/gw_analysis_tools/include/gwat/waveform_util.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gwat" TYPE FILE FILES "/opt/gw_analysis_tools/include/gwat/GWATConfig.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/opt/gw_analysis_tools/src/cmake_install.cmake")
  include("/opt/gw_analysis_tools/apps/cmake_install.cmake")
  include("/opt/gw_analysis_tools/data/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/opt/gw_analysis_tools/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
