
set (QTLMAP_ANALYSE_SOURCES
	 m_qtlmap_analyse_modlin_cox.f95  
	 m_qtlmap_analyse_multitrait.f95
	 m_qtlmap_analyse_multitrait_DA.f95
	 m_qtlmap_analyse_gen.f95  
	 m_qtlmap_analyse.F95       
	 m_qtlmap_analyse_lin_gen.f95   
	 m_qtlmap_analyse_modlin.f95 
	 m_qtlmap_analyse_modlin_ldla.f95
	 m_qtlmap_analyse_unitrait.f95
	 m_qtlmap_analyse_discret_unitrait.f95
	 m_qtlmap_incidence.F95
	 m_qtlmap_incidence_linear.F95
	 m_qtlmap_incidence_optim.F95
	 m_qtlmap_incidence_multi.f95
	 m_qtlmap_calcul_ic.f95
	 m_qtlmap_analyse_biallelic_farnir.f95
	)

IF (CUDA_IMPL)
find_package(CUDA REQUIRED)
set (CUDA_VERBOSE_BUILD on)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)

add_definitions(-DMANAGE_CUDA)
set (QTLMAP_ANALYSE_SOURCES_CUDA
     cuda_qtlmap_linear.cu
    )
    
set (QTLMAP_ANALYSE_SOURCES ${QTLMAP_ANALYSE_SOURCES} ${QTLMAP_ANALYSE_SOURCES_CUDA})

IF (NOT GENCODE_CUDA)
   set(GENCODE_CUDA "arch=compute_20,code=sm_20")
ENDIF (NOT GENCODE_CUDA)

set (CUDA_64_BIT_DEVICE_CODE true)
 
IF (TYPE_ARCH_CUDA MATCHES "_2")
set (CUDA_64_BIT_DEVICE_CODE true) 
ENDIF (TYPE_ARCH_CUDA MATCHES "_2")

SET (OPTIONS_CUDA -gencode=${GENCODE_CUDA})
# --compiler-bindir=${CMAKE_C_COMPILER}

IF (CUDA_BENCHMARK)
   SET (OPTIONS_CUDA ${OPTIONS_CUDA} -D_CUDA_HOST_TIME_PROF_=1)
ENDIF(CUDA_BENCHMARK)

IF (CUDA_PRINT_DEBUG)
   SET (OPTIONS_CUDA ${OPTIONS_CUDA} -D_CUDA_HOST_DEBUG_=1)
ENDIF (CUDA_PRINT_DEBUG)


CUDA_COMPILE(CUDA_QTLMAP_ANALYSE_SOURCES 
             ${QTLMAP_ANALYSE_SOURCES_CUDA}
             STATIC
             OPTIONS  ${OPTIONS_CUDA}
             RELEASE -O4 --use_fast_math
             DEBUG -g -G
             )

IF (CUDA_SP)
     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DCUDASP")
     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDASP")
     add_definitions(-DCUDASP)
ENDIF (CUDA_SP)

MESSAGE(STATUS " **************************************** ")
MESSAGE(STATUS " *             CUDA INFO                * ")
MESSAGE(STATUS " **************************************** ")

MESSAGE(STATUS " ** CUDA SOURCES :" ${CUDA_QTLMAP_ANALYSE_SOURCES})
IF (CUDA_SP)
  MESSAGE(STATUS " ** IMPL CUDA SIMPLE PRECISION ")
ELSE (CUDA_SP)
  MESSAGE(STATUS " ** IMPL CUDA DOUBLE PRECISION ")
ENDIF (CUDA_SP)

MESSAGE(STATUS " ** CUDA NVCC             : ${CUDA_NVCC}")
MESSAGE(STATUS " ** CUDA VERSION          : ${CUDA_VERSION_MAJOR}.${CUDA_VERSION_MINOR}")
MESSAGE(STATUS " ** CUDA LIBS             : ${CUDA_LIBRARIES}")
MESSAGE(STATUS " CODE GENERATION (GENCODE_CUDA) : ${GENCODE_CUDA} ")
MESSAGE(STATUS " **************************************** ")

set (QTLMAP_ANALYSE_SOURCES 
     ${QTLMAP_ANALYSE_SOURCES} 
     ${CUDA_QTLMAP_ANALYSE_SOURCES})
ENDIF (CUDA_IMPL)



include_directories(${QTLMAP_ALL_INCLUDE})
add_library (qtlmap_analyse ${QTLMAP_ANALYSE_SOURCES})      
add_dependencies(qtlmap_analyse qtlmap_data)
add_dependencies(qtlmap_analyse qtlmap_log)
add_dependencies(qtlmap_analyse qtlmap_output)
add_dependencies(qtlmap_analyse qtlmap_haplotype)

MESSAGE("** cmake analyse done. ")
