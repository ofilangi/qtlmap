set (QTLMAP_SOURCES
		 qtlmap.F95
	)
	 
set (TO_CARTHA_SOURCES
		 toCarthagene.f95
	)

set (QTLMAP_CHECK_SOURCES
		 qtlmap-check.f95
	)

include_directories (${QTLMAP_ALL_INCLUDE})
link_directories (${SRC_QTLMAP}/log)
link_directories (${SRC_QTLMAP}/tools)
link_directories (${SRC_QTLMAP}/data)
link_directories (${SRC_QTLMAP}/output/qtlmap_output)
link_directories (${SRC_QTLMAP}/input/qtlmap_input)
link_directories (${SRC_QTLMAP}/haplotype/qtlmap_haplotype)
link_directories (${SRC_QTLMAP}/analyse/qtlmap_analyse)


add_executable (${QTLMAP_PROG} ${QTLMAP_SOURCES})
target_link_libraries (${QTLMAP_PROG} qtlmap_log qtlmap_output qtlmap_input qtlmap_analyse qtlmap_haplotype qtlmap_tools qtlmap_data ${LIBS})

#si le main est dans .cpp
#set_target_properties(${QTLMAP_PROG} PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(${QTLMAP_PROG} stdc++)

#Le main est en fortran 
set_target_properties(${QTLMAP_PROG} PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS ${QTLMAP_PROG} RUNTIME DESTINATION "bin")

MESSAGE("** cmake progs done. ")
