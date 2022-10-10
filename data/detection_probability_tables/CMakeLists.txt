file(GLOB DET_PROB_TABLES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/data/detection_probability_tables/*.csv")
install(FILES ${DET_PROB_TABLES} DESTINATION "share/gwat/detection_probability_tables/")
