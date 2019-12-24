# Auteur : olivier.filangi@rennes.inra.fr
#
# Generating Release executable : 
# 
# cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release [./path-to-CMakeLists.txt] 
#
# Generating Release GPU C2050/C2075 executable : 
# cmake -DCUDA_IMPL=true -DGENCODE_CUDA="arch=compute_20,code=sm_20" [./path-to-CMakeLists.txt] 
#
# Directory to install QTLMAp   : cmake -DCMAKE_INSTALL_PREFIX=/path/to/binaries
# --------------------------- GPU OPTIONS ---------------------------------------------------------
# CUDA SDK ROOT DIR             : cmake -DCUDA_TOOLKIT_ROOT_DIR=/path
# CUDA ARCHITECTURE             : cmake -DGENCODE_CUDA="arch=compute_xx,code=sm_xx" xx=10,11,12,13,20
# SIMPLE PRECISION (CUDA)       : cmake -DCUDA_SP=true



# LAPACK/BLAS                   : cmake -DLAPACK_PATH=/path/lapack
# GMP                           : cmake -DGMP_PATH=/path/gmp
# ---------------------- DEVEL OPTIONS ---------------------------------------------------
# Generating Debug executable   : cmake -DCMAKE_BUILD_TYPE=Debug  
# use gprof                     : cmake -DGPROF=true
# view benchmark                : cmake -DBENCHMARK=true
# Force the compiler            : cmake -DCMAKE_Fortran_COMPILER=gfortran 
#                                       -DCMAKE_CXX_COMPILER=g++
#                                       -DCMAKE_C_COMPILER=gcc
# PRINT CUDA DEBUG INFO         : cmake -DCUDA_PRINT_DEBUG=true
# CUDA BENCHMARK                : cmake -DCUDA_BENCHMARK=true
# ---------------------- OLD OPTIONS -----------------------------------------------------
# Generating in NAG environment : cmake -DUSE_NAG_LIBRARY=true
# MPI                           : cmake -DMPI=true

# ---------------------- ARCH 32 --------------------------------------------------------------

# Compilation librairie NLOP seulement avec gcc comme compilateur C 
#                               : cmake -DCMAKE_ADDITIONAL_FLAGS="-m32"

#========================================================================================================================

#MPI
#Faire un setenv MPICH_F90 ifort / setenv MPICH_F90 gfortran / export MPICH_F90=ifort / export MPICH_F90=gfortran 
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/local/mpich/mpich2-1.1-gf44/bin/mpif90 -DMPI_INCLUDE_PATH=/local/mpich/mpich2-1.1-gf44 ../../
#cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Debug -DMPI_INCLUDE_PATH=/tools/cluster/openmpi-1.2.9/intel-11.0 ../../
#cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE=Debug ../../

# create a binary package       : make package
# create source package         : make package_source

#CMAKE REQUIREMENT
#----------------
cmake_minimum_required(VERSION 2.8 FATAL_ERROR) 

#DEBUG CMAKELISTS
#----------------
SET(CMAKE_VERBOSE_MAKEFILE false) 

#VERSIONNING
#-----------
SET(CPACK_PACKAGE_VERSION_MAJOR 1)
SET(CPACK_PACKAGE_VERSION_MINOR 0)
SET(CPACK_PACKAGE_VERSION_PATCH 0)

SET(QTLMAP_VERSION '${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}')
    
# CMake project file for QTLMap
#------------------------------
enable_language (Fortran)
PROJECT (QTLMap)

#base directory of include file for compilation step
set (INTERNAL_INCLUDE_QTLMAP ${QTLMap_SOURCE_DIR}/lib)

set (SRC_QTLMAP ${QTLMap_SOURCE_DIR}/src)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
	
SET(CPACK_PROJECT_NAME "QTLMap")
SET(CPACK_PACKAGE_NAME "qtlmap")
SET(CPACK_PACKAGE_FILE_NAME "qtlmap-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-bin-${CMAKE_SYSTEM_NAME}")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "QTLMap project")
SET(CPACK_PACKAGE_VENDOR "INRA - Animal Genetics Division")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README")
#SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/WELCOME")
		
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "QTLMap ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
IF(WIN32 AND NOT UNIX)
# There is a bug in NSI that does not handle full unix paths properly. Make
# sure there is at least one set of four (4) backlasshes.
#SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/\\\\InstallIcon.bmp")
SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\qtlmap.exe")
SET(CPACK_NSIS_DISPLAY_NAME "QTLMap") 

SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.inra.fr/qtlmap")
SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.inra.fr/qtlmap")
SET(CPACK_NSIS_CONTACT "olivier.filangi@rennes.inra.fr")
SET(CPACK_GENERATOR "NSIS")
SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
        set (CPACK_GENERATOR "TGZ")
		#SET(CPACK_STRIP_FILES "bin/qtlmap")
ENDIF(WIN32 AND NOT UNIX)

SET(CPACK_SOURCE_GENERATOR "TGZ")	
SET(CPACK_SOURCE_STRIP_FILES "") 
#filters...
SET(IGNORE_FILES "moc_" "\\\\.a$" "\\\\.moc" "\\\\.gmo" "ui_" "\\\\.gz" "\\\\.tar.Z" "\\\\.tar.bz2" "/doc/devel/" "\\\\.doc" "\\\\.odg" "\\\\.odt" "\\\\.tif" "\\\\.png" "gmon.out")
SET(IGNORE_FILES ${IGNORE_FILES} "\\\\.tar.xz" "\\\\~$" "\\\\kdevelop" "\\\\automoc" "Makefile$" "TODO" "CMakeCache.txt" )
SET(IGNORE_FILES ${IGNORE_FILES} "install_manifest.txt" "cmake_install.cmake" "CTest" "CPack" "/CVS/" "/CMakeTmp/" "\\\\.project" "\\\\.cproject")
SET(IGNORE_FILES ${IGNORE_FILES} "/CMakeFiles/;" "/_CPack_Packages/;" "/templates/;" "/\\\\.svn/;" "\\\\.swp$;" "\\\\.#;" "/#;" "/build/;" "/exemple/;" "/test/;" "/misc/;" "/hdoc/;")

#set filters
SET(CPACK_SOURCE_IGNORE_FILES ${IGNORE_FILES})
SET(CPACK_PACKAGE_EXECUTABLES "qtlmap" "qtlmap")

INCLUDE(CPack)
	
SET(NOW "`date +%d.%m.%Y-%H.%M.%S`")
add_definitions(-DQTLMAP_VERSION=${QTLMAP_VERSION})
#add_definitions(-DDATE_BUILD="'${NOW}'")

MESSAGE(STATUS "**************  Compilation on OS : ${CMAKE_SYSTEM_NAME}  *******************")
MESSAGE(STATUS "binary dir          : ${CMAKE_BINARY_DIR}")
MESSAGE(STATUS "source dir          : ${CMAKE_SOURCE_DIR}")
MESSAGE(STATUS "current binary dir  : ${CMAKE_CURRENT_BINARY_DIR}")
MESSAGE(STATUS "current source dir  : ${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE(STATUS "Using lib directory : ${INTERNAL_LIB_QTLMAP}")
MESSAGE(STATUS "OpenMP              : ${OPENMP}")
MESSAGE(STATUS "MPI                 : ${MPI_FOUND}")
   
include ("${CMAKE_SOURCE_DIR}/cmake_modules/CompileOptionsQTLMap.cmake")

set(RANDLIB_INCLUDE_DIR ${SRC_QTLMAP}/lib/randlib)
set(OPTIMIZATION_INCLUDE_DIR ${SRC_QTLMAP}/lib/optimization)
set(ORDERPACK_INCLUDE_DIR ${SRC_QTLMAP}/lib/orderpack-2.0)
set(SYMMAX2SAT_INCLUDE_DIR ${SRC_QTLMAP}/lib/symmax2sat)
set(MISCF77_INCLUDE_DIR ${SRC_QTLMAP}/lib/miscf77)
    
set (QTLMAP_DATA_INCLUDE ${SRC_QTLMAP}/data)
set (QTLMAP_LOG_INCLUDE ${SRC_QTLMAP}/log)
set (QTLMAP_TOOLS_INCLUDE ${SRC_QTLMAP}/tools)
set (QTLMAP_INPUT_INCLUDE ${SRC_QTLMAP}/input)
set (QTLMAP_OUTPUT_INCLUDE ${SRC_QTLMAP}/output)
set (QTLMAP_HAPLOTYPE_INCLUDE ${SRC_QTLMAP}/haplotype)
set (QTLMAP_ANALYSE_INCLUDE ${SRC_QTLMAP}/analyse)
	
set (ORDERPACK_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib/orderpack-2.0)
set (QTLMAP_DATA_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/data)
set (QTLMAP_LOG_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/log)
set (QTLMAP_TOOLS_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/tools)
set (QTLMAP_INPUT_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/input)
set (QTLMAP_OUTPUT_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/output)
set (QTLMAP_HAPLOTYPE_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/haplotype)
set (QTLMAP_ANALYSE_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/analyse)
	
set (QTLMAP_ALL_INCLUDE ${QTLMAP_ALL_INCLUDE} ${INTERNAL_LIBMOD_QTLMAP} ${QTLMAP_DATA_INCLUDE} ${QTLMAP_LOG_INCLUDE} ${QTLMAP_TOOLS_INCLUDE} ${QTLMAP_INPUT_INCLUDE} 
    ${QTLMAP_OUTPUT_INCLUDE} ${QTLMAP_HAPLOTYPE_INCLUDE} ${QTLMAP_ANALYSE_INCLUDE} ${VARYING_STRING_INCLUDE_DIR} ${ORDERPACK_INCLUDE_DIR} 
    ${ORDERPACK_MOD_DIR} ${QTLMAP_DATA_MOD_DIR} ${QTLMAP_LOG_MOD_DIR} ${QTLMAP_TOOLS_MOD_DIR} ${QTLMAP_INPUT_MOD_DIR} ${QTLMAP_OUTPUT_MOD_DIR}
    ${QTLMAP_HAPLOTYPE_MOD_DIR} ${QTLMAP_ANALYSE_MOD_DIR})

set(LIBS ${LIBS} orderpack randlib opti miscf77 miscf95 symmax2sat)

IF(NLOPT_ACTIVE)
  set(LIBS ${LIBS} nlopt)
ENDIF(NLOPT_ACTIVE)

IF (CUDA_IMPL)
find_package(CUDA REQUIRED)
set(LIBS ${LIBS} ${CUDA_LIBRARIES})
ENDIF(CUDA_IMPL)

FIND_PATH(GMP_INCLUDE_DIR gmp.h $ENV{HOME}/.local/include)

find_library(GMP_LIBRARY NAMES gmp PATHS ${GMP_PATH}/lib /usr/lib /usr/local/lib)

IF ( NOT GMP_LIBRARY )
              MESSAGE(FATAL_ERROR "$$$$$$$$ Can not find GMP library => install GMP : (ubuntu) sudo apt-get install libgmp3-dev / http://gmplib.org/")
ENDIF ( NOT GMP_LIBRARY )

IF (GMP_INCLUDE_DIR AND GMP_LIBRARY)
   SET(GMP_FOUND TRUE)
ENDIF (GMP_INCLUDE_DIR AND GMP_LIBRARY)


set(LIBS ${LIBS} ${GMP_LIBRARY})

# Recurse into the subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
     
     add_subdirectory (${QTLMap_SOURCE_DIR}/lib)
     
	 add_subdirectory (${SRC_QTLMAP}/log)
	 add_subdirectory (${SRC_QTLMAP}/tools)
	 add_subdirectory (${SRC_QTLMAP}/data)
	 add_subdirectory (${SRC_QTLMAP}/output)
	 add_subdirectory (${SRC_QTLMAP}/input)
	 add_subdirectory (${SRC_QTLMAP}/analyse)
	 add_subdirectory (${SRC_QTLMAP}/haplotype)
	 
	 set_target_properties(qtlmap_log PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_tools PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_data PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_output PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_input PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_analyse PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
     set_target_properties(qtlmap_haplotype PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
   
# build executables

    set (QTLMAP_PROG "qtlmap")
    set (QTLMAP_CHECK_PROG "qtlmap-check")
    set (EXECUTABLES ${CMAKE_CURRENT_BINARY_DIR}/src/${QTLMAP_PROG})
	add_subdirectory (${SRC_QTLMAP})
	set_target_properties(${QTLMAP_PROG} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGSF08}")
    
# default installation

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "Foo install prefix" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
MESSAGE(STATUS "Install directory                 : ${CMAKE_INSTALL_PREFIX}")

 
# install executables and scripts
install (FILES ${EXECUTABLES} DESTINATION bin)

#September 2013 = add unit Test with Fruit
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/testdev/unit)

MESSAGE("** cmake root done. ")
