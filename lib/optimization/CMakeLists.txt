set (OPTI_SOURCES
	 lbfgs-2.4.f
	) 
	 
add_library (opti ${OPTI_SOURCES})      


MESSAGE("** cmake optimization done. ")