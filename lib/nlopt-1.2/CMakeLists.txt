#add_subdirectory (api)
#d apres SHARED_VERSION_INFO
add_definitions(-DMAJOR_VERSION=2)
add_definitions(-DMINOR_VERSION=0)
add_definitions(-DBUGFIX_VERSION=2)

add_definitions(-UWITH_CXX)   

#  **************** STDINT.H --> TYPE UINT32 ***********************

find_file(STDINT NAMES stdint.h)

IF (STDINT-NOTFOUND)
    MESSAGE(STATUS "**** Can not found stdint.h ******")
    add_definitions(-UHAVE_STDINT_H)   
ELSE(STDINT-NOTFOUND)
    MESSAGE(STATUS ">> Found stdint.h : ${STDINT}")
    add_definitions(-DHAVE_STDINT_H)    
ENDIF(STDINT-NOTFOUND)

INCLUDE (CheckTypeSize) 
CHECK_TYPE_SIZE("unsigned int" SIZEOF_UNSIGNED_INT)
CHECK_TYPE_SIZE("unsigned long" SIZEOF_UNSIGNED_LONG)
add_definitions(-DSIZEOF_UNSIGNED_INT=${SIZEOF_UNSIGNED_INT})    
add_definitions(-DSIZEOF_UNSIGNED_LONG=${SIZEOF_UNSIGNED_LONG}) 
add_definitions(-DHAVE_UINT32_T) 
    
#SET(SIZEOF_UNSIGNED_INT ${CMAKE_SIZEOF_UNSIGNED_INT}) 
MESSAGE(STATUS "SIZEOF_UNSIGNED_INT: ${SIZEOF_UNSIGNED_INT}")
MESSAGE(STATUS "SIZEOF_UNSIGNED_LONG: ${SIZEOF_UNSIGNED_LONG}")


#  **************** TIME.H --> TYPE UINT32 ***********************
find_file(SYSTIME_H NAMES sys/time.h)

IF (SYSTIME)
    MESSAGE(FATAL-ERROR "**** Can not found sys/time.h ******")
    add_definitions(-UHAVE_SYS_TIME_H)   
ELSE(SYSTIME)
    MESSAGE(STATUS ">> Found sys/time.h : ${STDINT}")
    add_definitions(-DHAVE_SYS_TIME_H)      
ENDIF(SYSTIME)

find_file(TIME_H NAMES time.h)

IF (TIME_H-NOTFOUND)
    MESSAGE(FATAL-ERROR "**** Can not found time.h ******")
    add_definitions(-UTIME_WITH_SYS_TIME)   
ELSE(TIME_H-NOTFOUND)
    MESSAGE(STATUS ">> Found time.h : ${STDINT}")
    add_definitions(-DTIME_WITH_SYS_TIME)      
ENDIF(TIME_H-NOTFOUND)

# ******************************** ISNAN ************************
add_definitions(-DHAVE_ISNAN=1)

#C option
#-Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs -Wundef -Wconversion -Wmissing-prototypes -Wmissing-declarations"
#set (CMAKE_C_FLAGS "-Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -fopenmp  -pedantic -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs -Wundef -Wconversion -Wmissing-prototypes -Wmissing-declarations")

#set (CMAKE_C_FLAGS_RELEASE "-O2 ${CMAKE_C_FLAGS}")
#set (CMAKE_C_FLAGS_DEBUG "-g ${CMAKE_C_FLAGS}")
#CXX
#CXXFLAGS="$CXXFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wshadow -Wwrite-strings -Wredundant-decls -Wundef -Wconversion"
#set (CMAKE_CXX_FLAGS "-Wall -W -Wcast-qual -Wpointer-arith -fopenmp -Wcast-align -pedantic -Wshadow -Wwrite-strings -Wredundant-decls -Wundef -Wconversion")
#set (CMAKE_CXX_FLAGS_RELEASE "-O2 ${CMAKE_CXX_FLAGS}")
#set (CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")


set (NLOPT_SOURCES	qtlmap_interface.c)

aux_source_directory(api NLOPT_SOURCES)
aux_source_directory(auglag NLOPT_SOURCES)
aux_source_directory(bobyqa NLOPT_SOURCES)
aux_source_directory(cdirect NLOPT_SOURCES)
aux_source_directory(cobyla NLOPT_SOURCES)
aux_source_directory(crs NLOPT_SOURCES)
aux_source_directory(direct NLOPT_SOURCES)
aux_source_directory(isres NLOPT_SOURCES)
aux_source_directory(luksan NLOPT_SOURCES)
aux_source_directory(m4 NLOPT_SOURCES)
aux_source_directory(mlsl NLOPT_SOURCES)
aux_source_directory(mma NLOPT_SOURCES)
aux_source_directory(neldermead NLOPT_SOURCES)
aux_source_directory(newuoa NLOPT_SOURCES)
#aux_source_directory(octave NLOPT_SOURCES)
aux_source_directory(praxis NLOPT_SOURCES)
aux_source_directory(stogo NLOPT_SOURCES)
aux_source_directory(util NLOPT_SOURCES)


include_directories (api)
include_directories (auglag)
include_directories (bobyqa)
include_directories (cdirect)
include_directories (cobyla)
include_directories (crs)
include_directories (direct)
include_directories (isres)
#include_directories (lbfgs)
include_directories (luksan)
include_directories (m4)
include_directories (mlsl)
include_directories (mma)
include_directories (neldermead)
include_directories (newuoa)
include_directories (octave)
include_directories (praxis)
include_directories (stogo)
include_directories (util)

add_library (nlopt ${NLOPT_SOURCES})      

MESSAGE("** cmake nlopt done. ")
