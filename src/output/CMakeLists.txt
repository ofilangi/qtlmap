set (QTLMAP_OUTPUT_SOURCES
     m_qtlmap_output_handler.f95
     m_qtlmap_format_carthagene.f95		
	)

include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_output ${QTLMAP_OUTPUT_SOURCES})
add_dependencies(qtlmap_output qtlmap_data)
add_dependencies(qtlmap_output qtlmap_log)
add_dependencies(qtlmap_output qtlmap_tools)
add_dependencies(qtlmap_output qtlmap_haplotype)

MESSAGE("** cmake output done. ")

