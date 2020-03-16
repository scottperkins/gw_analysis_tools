
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/gwat/*.h")
file(GLOB SRC_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cpp")


add_library(OBJS OBJECT ${SRC_LIST})
set_property(TARGET OBJS PROPERTY  POSITION_INDEPENDENT_CODE ON)
set_property(TARGET OBJS PROPERTY CXX_STANDARD 11)
target_link_libraries( OBJS PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}" m)
target_include_directories( OBJS PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
target_include_directories( OBJS PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}" "${GSL_INCLUDE_DIR}")

add_library(gwat SHARED $<TARGET_OBJECTS:OBJS> ${HEADER_LIST})
add_library(gwat_static STATIC $<TARGET_OBJECTS:OBJS> ${HEADER_LIST})
set_property(TARGET gwat PROPERTY CXX_STANDARD 11)
set_property(TARGET gwat_static PROPERTY CXX_STANDARD 11)

#add_library(gwat SHARED ${SRC_LIST} ${HEADER_LIST})
#add_library(gwat_static STATIC ${SRC_LIST} ${HEADER_LIST})
target_include_directories( gwat PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
target_include_directories( gwat PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}")
target_link_libraries( gwat PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" "${GSL_LIBRARIES}" )

#target_link_libraries( gwat_static PUBLIC "${HDF5_LIBRARIES}" OpenMP::OpenMP_CXX "${ADOLC_LIBRARIES}" "${FFTW3_LIBRARIES}" )
#target_include_directories( gwat_static PUBLIC "${PROJECT_SOURCE_DIR}/include/gwat")
#target_include_directories( gwat_static PUBLIC "${OpenMP_INCLUDE_DIR}" "${ADOLC_INCLUDE_DIR}" "${HDF5_INCLUDE_DIR}" "${FFTW3_INCLUDE_DIR}")



#Build and name
set_target_properties(gwat_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set_target_properties(gwat PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set_target_properties(gwat_static PROPERTIES OUTPUT_NAME "gwat")

#Install
install(TARGETS gwat DESTINATION lib)
install(TARGETS gwat_static DESTINATION lib)
