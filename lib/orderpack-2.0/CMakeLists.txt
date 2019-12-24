set (ORDER_PACK_2_0_SOURCES
	ctrper.f90  indmed.f90  inspar.f90  mrgrnk.f90  refsor.f90  rnkpar.f90  unipar.f90  unista.f90  valnth.f90
	fndnth.f90  indnth.f90  inssor.f90  mrgref.f90  mulcnt.f90  refpar.f90  rinpar.f90  sort7.f90   uniinv.f90  unirnk.f90
	)
add_library (orderpack ${ORDER_PACK_2_0_SOURCES})

MESSAGE("** cmake orderpack done. ")
