set (QTLMAP_LOG_SOURCES
     m_qtlmap_log.F95	
	)

include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_log ${QTLMAP_LOG_SOURCES})      

MESSAGE("** cmake log done. ")