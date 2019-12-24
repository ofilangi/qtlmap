###############################################################################
#  Departement Genetique Animale
#  Author : Olivier.Filangi@rennes.inra.fr
#
#  Generic script to support several compiler for development
#
#  Option 
#    - -DOPENMP=true     : set true the preprocessor variable MANAGE_OMP in the source code 
#    - -DBENCHMARK=true  : set true the preprocessor variable BENCHMARK_VIEW in the source code
#    - -DCUDA_SP=true    : set true the preprocessor variable CUDA_SP in the source code
#    - -DCUDA_IMPL=true  : set true the preprocessor variable MANAGE_CUDA in the source code
#
#
###############################################################################

   
   set (MOD_QTLMAP ${CMAKE_CURRENT_BINARY_DIR}/mod)    
   
   #Ajout de la date pour la datation de la compilation de l executable
   set (ADD_OPTIONS "-DDATE_BUILD=\"'${NOW}'\"")
   
   #ajout de openmp 
   IF (OPENMP)
      add_definitions(-DMANAGE_OMP)
      set (ADD_OPTIONS "${ADD_OPTIONS} -DMANAGE_OMP")
   ENDIF (OPENMP)
   
   #Gestion du BENCHMARK de QTLMAP
   IF (BENCHMARK)
     set (ADD_OPTIONS "${ADD_OPTIONS} -DBENCHMARK_VIEW")
   ENDIF (BENCHMARK)
   
   #gestion des analyse cuda en simple precision, par defaut double precision 
   IF (CUDA_SP)
     add_definitions(-DCUDA_SP)
     set (ADD_OPTIONS "${ADD_OPTIONS} -DCUDA_SP")
   ENDIF (CUDA_SP)
   
   IF (CUDA_IMPL)
     add_definitions(-DMANAGE_CUDA)
     set (ADD_OPTIONS "${ADD_OPTIONS} -DMANAGE_CUDA")
   ENDIF (CUDA_IMPL)

   #ADDITIONAL FLAGS
   IF (NOT CMAKE_ADDITIONAL_FLAGS STREQUAL "")
     set (ADD_OPTIONS "${ADD_OPTIONS} ${CMAKE_ADDITIONAL_FLAGS}")
   ENDIF (NOT CMAKE_ADDITIONAL_FLAGS STREQUAL "")
   
# -------------------------------------------------------------------------  Manage CUBLAS librairy

   FIND_PATH(CUBLAS_CUBLASH_INCLUDE_DIR cublas.h ${CUBLAS_INCLUDE_PATH}/include /usr/local/cuda/include)
   FIND_PATH(CUBLAS_COMMONH_INCLUDE_DIR fortran_common.h ${CUBLAS_INCLUDE_PATH}/src /usr/local/cuda/src)

   IF ( CUBLAS_CUBLASH_INCLUDE_DIR AND CUBLAS_COMMONH_INCLUDE_DIR )
    SET(CUBLAS_FOUND TRUE)
   ENDIF ( CUBLAS_CUBLASH_INCLUDE_DIR AND CUBLAS_COMMONH_INCLUDE_DIR )

   IF (CUBLAS AND CUBLAS_FOUND)
     add_definitions(-DMANAGE_CUBLAS)
   ENDIF (CUBLAS AND CUBLAS_FOUND)

# ---------------------------------------------------------------------------
    
   MESSAGE(STATUS  "C Compiler: ${CMAKE_C_COMPILER}")
   if (CMAKE_C_COMPILER MATCHES "gcc") #gnu
    set(NLOPT_ACTIVE "true")
    set (ADD_OPTIONS "${ADD_OPTIONS} -DHAVE_NLOPT")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_OPTIONS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_OPTIONS}")
    
    MESSAGE(STATUS  " ** Active NLOPT ** ")
   else (CMAKE_C_COMPILER MATCHES "gcc")
    MESSAGE(STATUS  " ******************************* ")
    MESSAGE(STATUS  " ***** NLOPT IS NOT ACTIVE ***** ")
    MESSAGE(STATUS  " ******************************* ")
   endif (CMAKE_C_COMPILER MATCHES "gcc")
       
   IF (WIN32)
    add_definitions(-DWINDOWS)
    set (ADD_OPTIONS "${ADD_OPTIONS} -DWINDOWS")
   ELSE (WIN32)
    add_definitions(-DLINUX)
    set (ADD_OPTIONS "${ADD_OPTIONS} -DLINUX")
   ENDIF (WIN32)
   
   get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
   MESSAGE(STATUS  "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
      
# make sure that the default is a RELEASE
	if (NOT CMAKE_BUILD_TYPE)
	  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
	endif (NOT CMAKE_BUILD_TYPE)
    
    
   #Gestion MPI
   #-----------
   IF (Fortran_COMPILER_NAME MATCHES "mpif90")
       MESSAGE(STATUS  "MPI Found")
       FIND_PATH(MPI_INCLUDE_DIR mpi.h ${MPI_INCLUDE_PATH}/include /usr/include/mpi /usr/local/include/mpi)
      # FIND_LIBRARY(MPI_LIBRARY NAMES mpi PATH  ${MPI_INCLUDE_PATH}/lib /usr/lib /usr/local/) 

       IF (MPI_INCLUDE_DIR)
          SET(MPI_FOUND TRUE)
       ENDIF (MPI_INCLUDE_DIR)


       IF (MPI_FOUND)
         IF (NOT MPI_FIND_QUIETLY)
          MESSAGE(STATUS "Found MPI: ${MPI_LIBRARY}")
         ENDIF (NOT MPI_FIND_QUIETLY)
       ELSE (MPI_FOUND)
         IF (MPI_FIND_REQUIRED)
       #   MESSAGE(FATAL_ERROR "Could not find MPI")
         ENDIF (MPI_FIND_REQUIRED)
       ENDIF (MPI_FOUND)

      # include_directories(${MPI_INCLUDE_PATH})
       set (QTLMAP_ALL_INCLUDE "${QTLMAP_ALL_INCLUDE} ${MPI_INCLUDE_DIR}" )
       #set(LIBS ${LIBS} ${MPI_LIBRARY})
       MESSAGE(STATUS "==>INCLUDE MPI:${MPI_INCLUDE_DIR}")
       #MESSAGE(STATUS "==>LIB MPI:${MPI_LIBRARY}")

       link_directories(${MPI_LIBRARY})
     #  set(CMAKE_Fortran_COMPILER "mpif90")
       set(Fortran_COMPILER_NAME $ENV{MPICH_F90})
       
       IF (DEFINED Fortran_COMPILER_NAME) 
         MESSAGE(STATUS  "*** Compilo MPI :${Fortran_COMPILER_NAME} *****")
       ELSE ( DEFINED Fortran_COMPILER_NAME)
         MESSAGE(FATAL_ERROR  "**** Set environment variable MPICH_F90 (setenv MPICH_F90 <compiler>) | export MPICH_F90=<compiler>) *** ")
       ENDIF (DEFINED Fortran_COMPILER_NAME)

     #  SET(ENV{MPICH_F90} "${Fortran_COMPILER_NAME}")
     #  SET(ENV{MPICH_F90LINKER} "${Fortran_COMPILER_NAME}")
       add_definitions(-DMPI_ACTIVE)
       set (ADD_OPTIONS "${ADD_OPTIONS} -DMPI_ACTIVE")
 #option pour debugging...a mettre en place
 #-mpe=mpicheck
 #-mpe=mpitrace
 #    ELSE (MPI_FOUND)
 #        MESSAGE(STATUS  "MPI Not Found (required to compile modules)")
 #    ENDIF (MPI_FOUND)
    ENDIF (Fortran_COMPILER_NAME MATCHES "mpif90")
    
    
#----------------------- NAG ---------------------------------------------

IF (USE_NAG_LIBRARY)
       find_library(NAG_LIBRARY NAMES nag_nag nag nag64)
       IF ( NOT NAG_LIBRARY )
              MESSAGE(FATAL_ERROR "$$$$$$$$ Can not find NAG library $$$$$$$$$ ")
       ELSE ( NOT NAG_LIBRARY )
              MESSAGE(STATUS "******************************************* ")
              MESSAGE(STATUS "************** NAG library **************** ")
              MESSAGE(STATUS "******************************************* ")
              MESSAGE(STATUS "Find NAG: ${NAG_LIBRARY}")
              #add_library(nag_nag SHARED ${NAG_LIBRARY})
              find_library(PHREAD_LIBRARY NAMES pthread)
              set(LIBS ${LIBS} ${NAG_LIBRARY} ${PHREAD_LIBRARY})
              IF (NOT PHREAD_LIBRARY)
                MESSAGE (FATAL_ERROR "Can not find library pthread")
              ENDIF (NOT PHREAD_LIBRARY)
              add_definitions(-DHAVE_LIBNAG)
       ENDIF ( NOT NAG_LIBRARY )
ENDIF (USE_NAG_LIBRARY)
    
#----------------------- GFORTRAN ----------------------------------------	
	if (Fortran_COMPILER_NAME MATCHES "gfortran" OR Fortran_COMPILER_NAME MATCHES "f95") #gnu
	   execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion OUTPUT_VARIABLE GFORTRAN_VERSION)
	   string(REGEX MATCHALL "[0-9]+" GFORTRAN_VERSION_COMPONENTS ${GFORTRAN_VERSION})
           list(GET GFORTRAN_VERSION_COMPONENTS 0 GFORTRAN_MAJOR)

           message(STATUS "Gfortran version :  ${GFORTRAN_MAJOR}")
	   #preprocessor
           set (ADD_OPTIONS "${ADD_OPTIONS} -cpp")

           # DEPENDANCES : LAPACK / BLAS
           # 20/06/2011 
           # 
           find_library(LAPACK_LIBRARY NAMES lapack PATHS ${LAPACK_PATH}/lib /usr/lib /usr/local/lib)
           find_library(BLAS_LIBRARY NAMES blas PATHS ${LAPACK_PATH}/lib /usr/lib /usr/local/lib)

           IF ( NOT LAPACK_LIBRARY )
              MESSAGE(FATAL_ERROR "$$$$$$$$ Can not find LAPACK library $$$$$$$$$ => install lapack : http://www.netlib.org/lapack/")
           ENDIF ( NOT LAPACK_LIBRARY )

           IF ( NOT BLAS_LIBRARY )
              MESSAGE(FATAL_ERROR "$$$$$$$$ Can not find Blas library $$$$$$$$$ => install lapack : http://www.netlib.org/lapack/")
           ENDIF ( NOT BLAS_LIBRARY )

           set(LIBS ${LIBS} ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

	   #IF ( GFORTRAN_MAJOR EQUAL 4 ) #only major version 4
	     IF (GPROF)
           set (ADD_OPTIONS "${ADD_OPTIONS} -pg")
         ENDIF (GPROF)
       #
	        set (CMAKE_Fortran_MODULE_DIRECTORY ${MOD_QTLMAP} )
	        #set (ADD_OPTIONS "${ADD_OPTIONS} -DGNU_COMP -DMANAGE_OMP -fopenmp")
          set (ADD_OPTIONS "${ADD_OPTIONS} -DGNU_COMP")
   	      set (CMAKE_Fortran_FLAGS_RELEASE "${ADD_OPTIONS} -O3 -fexternal-blas")
          set (CMAKE_Fortran_FLAGS_DEBUG   "${ADD_OPTIONS} -g3 -O0 -ggdb -fbacktrace -fdump-core --coverage -fstack-usage -fcheck=bounds")
          set (CMAKE_Fortran_FLAGSF77  "")
          set (CMAKE_Fortran_FLAGSF90  "")
          set (CMAKE_Fortran_FLAGSF95  "-std=f95")
          set (CMAKE_Fortran_FLAGSF08  "-std=f2008")
	      
	   #ENDIF ( GFORTRAN_MAJOR EQUAL 4 )
#----------------------- IFORT ----------------------------------------
	elseif (Fortran_COMPILER_NAME MATCHES "ifort") # intel
    
        find_package(MKL REQUIRED)
        IF(MKL_FOUND)
        include_directories(${MKL_INCLUDE_DIRS})
        #link_directories(${MKL_LIBRARIES})
        SET(LIBS ${LIBS} ${MKL_LIBRARIES})
   ELSE()
    message(FATAL_ERROR "--> can not find MKL library")
   ENDIF()

       set (CMAKE_Fortran_MODULE_DIRECTORY ${MOD_QTLMAP} )
       FIND_PACKAGE(IntelQTLMapPref) 
           set (ADD_OPTIONS "${ADD_OPTIONS} -DINTEL_COMP -openmp -openmp-lib=compat -DMANAGE_OMP")
	   set (CMAKE_Fortran_FLAGS_RELEASE "${ADD_OPTIONS} -O3 -heap-array")
           set (CMAKE_Fortran_FLAGS_DEBUG   "${ADD_OPTIONS} -O0 -g -fp-stack-check -traceback -debug all")	    

	   set (CMAKE_Fortran_FLAGSF77 "-f77rtl")
           set (CMAKE_Fortran_FLAGSF90 "-fpp -stand f90 -free")
	   set (CMAKE_Fortran_FLAGSF95 "-fpp -stand f95 -free")
           set (CMAKE_Fortran_FLAGSF08 "-fpp -stand f08 -free")
	   
#----------------------- G95 ----------------------------------------
	elseif (Fortran_COMPILER_NAME MATCHES "g95") 
	    
	   set (CMAKE_Fortran_FLAGS_RELEASE "${ADD_OPTIONS} -O3  -fno-second-underscore  ")
	   set (CMAKE_Fortran_FLAGS_DEBUG   "${ADD_OPTIONS} -O0  -fno-second-underscore -g -fbounds-check -ftrace=frame -fzero")
	   set (CMAKE_Fortran_FLAGSF77 "")
           set (CMAKE_Fortran_FLAGSF90 "")
	   set (CMAKE_Fortran_FLAGSF95 "-cpp")
	   set (CMAKE_Fortran_FLAGSF08 "-cpp")

	    #linker option
       set(LIBS ${LIBS} -lf95) 
       set(LIBS ${LIBS} -lgfortran)
# ------------------------ ERROR : NO COMPILER FOUND ------------------
    else (Fortran_COMPILER_NAME MATCHES "gfortran" OR Fortran_COMPILER_NAME MATCHES "f95")
      message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
      message ("CMAKE_CXX_COMPILER full path: " ${CMAKE_CXX_COMPILER})
	  message ("Cxx compiler: " ${CXX_COMPILER_NAME})
	  message (FATAL_ERROR "The compiler no matches with this configuration file.")
    endif (Fortran_COMPILER_NAME MATCHES "gfortran" OR Fortran_COMPILER_NAME MATCHES "f95")
    

 #affichage
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	  MESSAGE(STATUS "Build type : Debug ")
	  MESSAGE(STATUS ${CMAKE_Fortran_FLAGS_DEBUG})
	else (CMAKE_BUILD_TYPE STREQUAL "Debug")
	  MESSAGE(STATUS "Build type : Release ")
	  MESSAGE(STATUS "Compiler Options :${CMAKE_Fortran_FLAGS_RELEASE}")
	endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
	
	MESSAGE(STATUS "F77 : ${CMAKE_Fortran_FLAGSF77}")
        MESSAGE(STATUS "F90 : ${CMAKE_Fortran_FLAGSF90}")
	MESSAGE(STATUS "F95 : ${CMAKE_Fortran_FLAGSF95}")
        MESSAGE(STATUS "F08 : ${CMAKE_Fortran_FLAGSF08}")
	
	MESSAGE("** general compiler configuration done.")
        MESSAGE("** Librairies = ${LIBS}")
	
	
	
    
