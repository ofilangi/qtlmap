SET(FRUIT_SRC 
	fruit/fruit_util.f90
    fruit/fruit.f90)

add_library (fruit ${FRUIT_SRC}) 
set_target_properties(fruit PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF90}")


set (TEST_SRC
	        main.f95
	        src/data/test_type_map.f95
	        src/data/test_type_lrt_solution.f95
)

add_executable(testdata ${TEST_SRC})
set_target_properties(testdata PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(testdata PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
target_link_libraries (testdata qtlmap_data qtlmap_log ${LIBS} fruit)