cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_FC_COMPILER mpif90)

project(abcd CXX C Fortran)

# Dependencies
set(EIGEN_DIR /usr/local/include/eigen3)
set(PETSC_DIR /home/knuthy/clones/petsc-3.3-p2/)
set(MPI_DIR /opt/mpich2)
set(BLAS_DIR /opt/OpenBLAS)

set(
  PETSC_INCLUDE_DIR 
  ${PETSC_DIR}/include 
  ${PETSC_DIR}/linux/include
  )
  
set(PETSC_LIBRARY_DIR ${PETSC_DIR}/linux/lib)

set(MPI_INCLUDE_DIR ${MPI_DIR}/include)
set(MPI_LIBRARY_DIR ${MPI_DIR}/lib)

# the library
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE}/lib)


link_directories(${Boost_LIBRARY_DIRS} ${PETSC_LIBRARY_DIR} ${MPI_LIBRARY_DIR} ${BLAS_DIR})
include_directories(include ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_DIR} ${PETSC_INCLUDE_DIR} ${EIGEN_DIR})

file(
	GLOB_RECURSE
	source_files
	src/*.f
	src/*
	include/*
)

add_library(abcd STATIC ${source_files})


add_subdirectory(test)
