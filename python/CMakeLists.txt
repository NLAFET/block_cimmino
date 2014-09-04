cmake_minimum_required(VERSION 2.6)

set(BUILD_TYPE "Release")

include(../abcdCmake.in)
project(abcdpy)

set(CMAKE_CXX_FLAGS "-fopenmp")
set(CMAKE_C_FLAGS "-fopenmp")

#convert ${ABCD} to absolute
get_filename_component(ABCD_DIR "${ABCD}" ABSOLUTE)

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

add_subdirectory("${CMAKE_SOURCE_DIR}/numpy")

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
  ${ABCD_DIR}/lib
  ${ABCD_DIR}/lib/sparselib/lib/
  ${ABCD_DIR}/lib/mumps/lib/
  ${ABCD_DIR}/build/lib
  ${Boost_LIBRARY_DIRS}
  ${MUMPS_LIBRARY_DIR}
  ${BLAS_LAPACK_SCALAPACK_DIRS}
  ${MPI_LIB_DIR}
  ${LIB_DIRS}
  )

message(${ABCD_DIR})

include_directories(
  ${ABCD_DIR}/include
  ${ABCD_DIR}/lib/sparselib/include
  ${ABCD_DIR}/lib/sparselib/mv/include
  ${ABCD_DIR}/lib/mumps/include
  ${Boost_INCLUDE_DIRS}
  ${BOOST_ROOT}
  ${MUMPS_INCLUDE_DIR}
  ${MPI_INC_DIR}
  ${CMAKE_SOURCE_DIR}/numpy/
  )

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(abcdpy SHARED
  python_interface.cpp
)
set_target_properties(abcdpy PROPERTIES PREFIX "")

target_link_libraries(abcdpy
  abcd
  boost_numpy
  ${LIBS}
  )
