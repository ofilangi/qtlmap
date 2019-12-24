
if (MKL_LIBRARIES)
  set(MKL_FIND_QUIETLY TRUE)
endif (MKL_LIBRARIES)

if(CMAKE_MINOR_VERSION GREATER 4)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

find_library(MKL_CORE
  mkl_core
  HINTS
  $ENV{MKLLIB}
  /opt/mkl-*/lib/em64t
  /opt/intel/mkl/*/lib/em64t
  /opt/intel/mkl/lib/intel64/
  /opt/intel/Compiler/*/*/mkl/lib/em64t
  ${LIB_INSTALL_DIR}
)

if(MKL_CORE)
  set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_CORE})
endif()

find_library(MKL_INTEL_LP64
  mkl_intel_lp64
  HINTS
  $ENV{MKLLIB}
  /opt/mkl-*/lib/em64t
  /opt/intel/mkl/*/lib/em64t
  /opt/intel/mkl/lib/intel64/
  /opt/intel/Compiler/*/*/mkl/lib/em64t
  ${LIB_INSTALL_DIR}
)

if(MKL_INTEL_LP64)
  set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_INTEL_LP64})
endif()

find_library(MKL_SEQUENTIAL
  mkl_sequential
  HINTS
  $ENV{MKLLIB}
  /opt/mkl-*/lib/em64t
  /opt/intel/mkl/*/lib/em64t
  /opt/intel/mkl/lib/intel64/
  /opt/intel/Compiler/*/*/mkl/lib/em64t
  ${LIB_INSTALL_DIR}
)


if(MKL_SEQUENTIAL)
  set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_SEQUENTIAL})
endif()


endif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

find_path(MKL_INCLUDE_DIR
  NAMES mkl.h
  HINTS
  $ENV{MKLINCLUDE}
  /opt/mkl-*/
  /opt/intel/mkl/*/
  /opt/intel/Compiler/*/*/mkl/
  ${CMAKE_INSTALL_PREFIX}/include
  PATH_SUFFIXES include
  )

endif(CMAKE_MINOR_VERSION GREATER 4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIR)

mark_as_advanced(MKL_LIBRARIES MKL_INCLUDE_DIR)
