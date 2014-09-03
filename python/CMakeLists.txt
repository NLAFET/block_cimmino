cmake_minimum_required(VERSION 2.6)

set(BUILD_TYPE "Debug")

include(../abcdCmake.in)
project(abcdpy)

set(CMAKE_CXX_FLAGS "-fopenmp")
set(CMAKE_C_FLAGS "-fopenmp")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package (Boost COMPONENTS mpi serialization python REQUIRED)
find_package (PythonLibs REQUIRED)

if(PYTHONLIBS_FOUND)
  include_directories("${PYTHON_INCLUDE_DIRS}")
else()
  message("PythonLibs didn't work")
endif()

# the library
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

if(NOT PATOH MATCHES OFF)
  set(LIBS ${LIBS} patoh)
endif()

set(LIBS ${LIBS}
  dmumps
  mumps_common

  ${Boost_LIBRARIES}
  ${BLAS_LAPACK_SCALAPACK_LIBS}
  ${MPI_LIBRARIES}

  sparse
  spblas
  mv
  gfortran
  )

link_directories(
  ${ABCD}/lib
  ${ABCD}/lib/sparselib/lib/
  ${ABCD}/lib/mumps/lib/
  ${ABCD}/build/${BUILD_TYPE}/lib
  ${Boost_LIBRARY_DIRS}
  ${MUMPS_LIBRARY_DIR}
  ${BLAS_LAPACK_SCALAPACK_DIRS}
  ${MPI_LIB_DIR}
  ${LIB_DIRS}
  )

include_directories(
  ${ABCD}/include
  ${ABCD}/lib/sparselib/include
  ${ABCD}/lib/sparselib/mv/include
  ${ABCD}/lib/mumps/include
  ${Boost_INCLUDE_DIRS}
  ${BOOST_ROOT}
  ${MUMPS_INCLUDE_DIR}
  ${MPI_INC_DIR}
  )


add_library(abcdpy SHARED
  python_interface.cpp
)
set_target_properties(abcdpy PROPERTIES PREFIX "")

target_link_libraries(abcdpy
  abcd_shared
  ${LIBS}
  )
