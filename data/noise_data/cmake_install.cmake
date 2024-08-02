# Install script for directory: /opt/gw_analysis_tools/data/noise_data

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gwat/noise_data" TYPE FILE FILES
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/AdLIGODesign.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/AdLIGOMidHigh.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/AplusDesign.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/AplusDesign_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/CE1_strain.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/CE1_strain_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/CE2_strain.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/CE2_strain_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/ET-0000A-18_ETDSensitivityCurveTxtFile.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/ETDXylophoneDwyer.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/Hanford_O2_Strain.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/Voyager.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/aligo_O4high.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/aligo_O4high_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O4high_NEW.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O4high_NEW_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O5high_NEW.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O5high_NEW_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O5low_NEW.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/avirgo_O5low_NEW_smoothed.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/kagra_128Mpc.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/kagra_25Mpc.csv"
    "/opt/gw_analysis_tools/data/noise_data/currently_supported/kagra_80Mpc.csv"
    )
endif()

