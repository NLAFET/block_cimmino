cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

#Configuration du projet
project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})
#set(MUMPS_DIR /opt/MUMPS_4.10.0)
#set(MPI_DIR /usr/lib/openmpi)
#set(BLAS /opt/OpenBLAS)
#set(LAPACK /opt/lapack)
#set(SLAPACK /opt/scalapack-2.0.1)
#set(BLACS /opt/BLACS/LIB)
#set(METIS /opt/ParMetis-3.1.1)
#set(SCOTCH /opt/scotch_5.1.11_esmumps)

find_package (Boost COMPONENTS mpi REQUIRED)

link_directories(
  ../${CMAKE_BUILD_TYPE}
  ${Boost_LIBRARY_DIRS}
  /usr/lib/openblas-base
#  ${MUMPS_DIR}/lib
#  ${MPI_DIR}
#  ${BLAS}
#  ${LAPACK}
#  ${SLAPACK}
#  ${BLACS}
#  ${METIS}
#  ${SCOTCH}/lib
  )
include_directories(../include ${Boost_INCLUDE_DIRS})

add_executable(
	test_abcd
	src/main.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
	test_abcd
  abcd 
	${Boost_LIBRARIES}
	
  lapack
  
  openblas
  
  scalapack-openmpi
  
  blacsF77init-openmpi
  blacsCinit-openmpi
  blacs-openmpi
  blacsF77init-openmpi
  blacsCinit-openmpi
  
  parmetis
  
  metis
  
  esmumps
  scotch
  scotcherr
  
	dmumps
	mumps_common
	
	mpi_f77
  mpi_f90
)
