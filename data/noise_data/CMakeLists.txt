file(GLOB NOISE_DATA CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/data/noise_data/currently_supported/*.csv")
install(FILES ${NOISE_DATA} DESTINATION "share/gwat/noise_data/")
