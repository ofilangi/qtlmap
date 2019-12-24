set (QTLMAP_HAPLOTYPE_SOURCES
     m_qtlmap_haplotype.f95  
	 m_qtlmap_haplotype_util.f95  
	 m_qtlmap_haplotype_V1.f95  
	 m_qtlmap_haplotype_V2.f95
	 m_qtlmap_haplotype_ldla.f95 
	 m_qtlmap_haplotype_external.f95 
     m_qtlmap_phase_offspring.f95 
	)

include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_haplotype ${QTLMAP_HAPLOTYPE_SOURCES})      
add_dependencies(qtlmap_haplotype qtlmap_data)
add_dependencies(qtlmap_haplotype qtlmap_log)
add_dependencies(qtlmap_haplotype qtlmap_tools)

MESSAGE("** cmake haplotype done. ")
