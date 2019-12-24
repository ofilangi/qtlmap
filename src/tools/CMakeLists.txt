set (QTLMAP_TOOLS_SOURCES
	m_qtlmap_math.F95
	m_qtlmap_solve_linear.F95
	m_qtlmap_optimization.F95
	m_qtlmap_isymmax2sat.f95   
)

include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_tools ${QTLMAP_TOOLS_SOURCES})
add_dependencies(qtlmap_tools qtlmap_data)     
add_dependencies(qtlmap_tools qtlmap_log)
ADD_DEPENDENCIES(qtlmap_tools ${LIBS})

MESSAGE("** cmake tools done. ")
