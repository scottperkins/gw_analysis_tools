
add_executable(util "src/test_util.cpp" )
target_link_libraries(util PUBLIC gwat)
target_link_libraries( util PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}")
target_include_directories( util PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
target_include_directories( util PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}" "${GSL_INCLUDE_DIR}")

add_executable(fishers "src/test_fishers.cpp" )
target_link_libraries(fishers PUBLIC gwat)
target_link_libraries( fishers PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}")
target_include_directories( fishers PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
target_include_directories( fishers PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}" "${GSL_INCLUDE_DIR}")

add_executable(snr "src/test_SNR.cpp" )
target_link_libraries(snr PUBLIC gwat)
target_link_libraries( snr PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}")
target_include_directories( snr PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
target_include_directories( snr PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}" "${GSL_INCLUDE_DIR}")

#Waveform testing needs LALsuite, which needs more tweaking
#add_executable(waveform "src/test_waveform.cpp" )
#target_link_libraries(waveform PUBLIC gwat)
#target_link_libraries( waveform PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}")
#target_include_directories( waveform PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
#target_include_directories( waveform PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}" "${GSL_INCLUDE_DIR}")
