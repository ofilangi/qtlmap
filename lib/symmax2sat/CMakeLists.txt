set (SYMMAX_SRC	tb2abstractconstr.cpp  tb2constraint.cpp   tb2system.cpp
tb2alldiffconstr.cpp   tb2domain.cpp                   tb2naryconstr.cpp     tb2ternaryconstr.cpp
tb2arithmetic.cpp      tb2enumvar.cpp                  tb2pedigree.cpp       tb2vac.cpp
tb2bep.cpp             tb2flowbasedconstr.cpp          tb2queue.cpp          tb2vacutils.cpp
   tb2globalcardinalityconstr.cpp  tb2randomgen.cpp      tb2variable.cpp
tb2boostgraph.cpp      tb2globalconstr.cpp             tb2reader.cpp         tb2wcsp.cpp
tb2btd.cpp             tb2graph.cpp                    tb2regularconstr.cpp
tb2btqueue.cpp         tb2haplotype.cpp                tb2sameconstr.cpp
tb2clusters.cpp        tb2intervar.cpp                 tb2solver.cpp tb2binconstr.cpp    )
	 
	 
#-fleading-underscore
add_definitions(-DWCSPFORMATONLY)
add_definitions(-DNARYCHAR)
add_definitions(-DNDEBUG)
add_definitions(-DLONGLONG_COST)
add_definitions(-DLINUX)
add_definitions(-DWIDE_STRING)
add_definitions(-DLONGDOUBLE_PROB )

# svn checkout svn://scm.mulcyber.toulouse.inra.fr/svnroot/toulbar2/trunk/toulbar2 
#modifcation of toulbar cp {TOULBAR}/src/* {QTLMAP}/symmax2sat/

# Modif in tb2solver.hpp et add in tb2solver.cpp
# ----------------------------------------------------------------------------------------------
#  int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol)
#
# ----------------------------------------------------------------------------------------------
# Atention int n, int m ===> int *n, int *m,
#
#  extern "C" int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost, int *sol);
# ----------------------------------------------------------------------------------------------


#int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost, int *sol)
#{
#  // select verbosity during search
#  ToulBar2::verbose = 0;
#  ToulBar2::uai=1;
#  //  ToulBar2::elimDegree = -1;
#
#  initCosts(MAX_COST);
#  Solver solver(STORE_SIZE, MAX_COST);
#
#  ToulBar2::startCpuTime = cpuTime();
#
#  return solver.solve_symmax2sat(*n , *m, posx, posy, cost, sol);
#}

FIND_PATH(GMP_PATH_DIR gmp.h ${GMP_PATH}/include /usr/include /usr/local/include)

IF ( NOT GMP_PATH_DIR )
     MESSAGE (FATAL_ERROR "Can not find GMP library  cmake -DGMP_PATH=/path/gmp" )  
ENDIF (NOT GMP_PATH_DIR )

set(CMAKE_CXX_FLAGS "-Warray-bounds -static-libgcc -fno-leading-underscore -static -Wall -I${GMP_PATH_DIR}")

add_library (symmax2sat ${SYMMAX_SRC})      


MESSAGE("** cmake symmax2sat done. ")
