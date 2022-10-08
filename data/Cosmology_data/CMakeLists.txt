file(GLOB COSMOLOGY_DATA CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/data/Cosmology_data/*.csv")
install(FILES ${COSMOLOGY_DATA} DESTINATION "share/gwat/Cosmology_data/")
