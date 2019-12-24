set (OPTI_SOURCES_F77
	 plott.F
     prob.f
	)

set (OPTI_SOURCES_F95
     dcdflib.f95
	)

add_library (miscf77 ${OPTI_SOURCES_F77})      

add_library (miscf95 ${OPTI_SOURCES_F95})      

MESSAGE("** cmake misc done. ")
