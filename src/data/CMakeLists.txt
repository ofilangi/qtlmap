set (QTLMAP_DATA_SOURCES
     m_qtlmap_constant.f95
     m_qtlmap_base.F95
     m_qtlmap_type_genealogy.f95
     m_qtlmap_type_map.f95
     m_qtlmap_type_phenotype.f95
     m_qtlmap_type_genotype.f95
     m_qtlmap_type_lrt_solution.f95
     m_qtlmap_type_incidence_solution.f95    
     m_qtlmap_type_dataset.f95
     m_qtlmap_type_dataset_trio.f95
     m_qtlmap_type_pdd.f95
     m_qtlmap_type_cli.F95
     m_qtlmap_type_parameter.f95
     m_qtlmap_type_parameter_simulation.f95
     m_qtlmap_types.f95
	)


include_directories(${QTLMAP_ALL_INCLUDE})
#SET_SOURCE_FILES_PROPERTIES(m_qtlmap_types.f95 PROPERTIES OBJECT_DEPENDS m_qtlmap_type_map.f95.o)
add_library (qtlmap_data ${QTLMAP_DATA_SOURCES}) 
add_dependencies(qtlmap_data qtlmap_log)
MESSAGE("** cmake data done. ")



