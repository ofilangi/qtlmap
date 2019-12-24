set (QTLMAP_INPUT_SOURCES
     m_qtlmap_simulation.f95  
     m_qtlmap_phenotype.f95
	 m_qtlmap_genealogy.f95  
	 m_qtlmap_genotype.f95  
	 m_qtlmap_map.f95
	 m_qtlmap_factory_dataset.f95
	)

include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_input ${QTLMAP_INPUT_SOURCES})
add_dependencies(qtlmap_input qtlmap_data)
add_dependencies(qtlmap_input qtlmap_output)
add_dependencies(qtlmap_input randlib)

MESSAGE("** cmake input done. ")

