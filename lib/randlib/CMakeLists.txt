set (RANDLIB_SOURCES
	 getcgn.f inrgcm.f  mltmod.f ranf.f setgmn.f ignlgi.f lennob.f  phrtsd.f sdot.f snorm.f genmn.f initgn.f qrgnin.f setall.f spofa.f
	)

add_library (randlib ${RANDLIB_SOURCES})      


MESSAGE("** cmake randlib done. ")