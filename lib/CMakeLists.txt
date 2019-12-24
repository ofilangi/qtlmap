   add_subdirectory (miscf77)
   add_subdirectory (randlib)
   add_subdirectory (optimization)
   add_subdirectory (orderpack-2.0)
   add_subdirectory (symmax2sat)
   
#la fonction d optimisation est une fonction definit dans une fonction...ce n est pas C ansi, seul gcc accept cette declaration
if (CMAKE_C_COMPILER MATCHES "gcc") #gnu
   add_subdirectory (nlopt-1.2)
   add_definitions(-DHAVE_NLOPT)
endif (CMAKE_C_COMPILER MATCHES "gcc")
 
   set_target_properties(miscf77 PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF77}")
   set_target_properties(miscf95 PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF90}")
   set_target_properties(randlib PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF77}")
   set_target_properties(opti PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF77}")
   set_target_properties(orderpack PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF95}")

MESSAGE("** cmake libs done. ")