cmake_minimum_required(VERSION 2.6)

set(BUILD_TYPE "Release")

include(../abcdCmake.in)
project(example)

find_package (Boost COMPONENTS mpi serialization REQUIRED)

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
  ${ABCD}/build/Release/lib
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


add_executable(
  example
  example.cpp
  )

###Configuration de l'édition de liens
target_link_libraries(
  example
  abcd 
  ${LIBS}
  )
