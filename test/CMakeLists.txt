cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpic++)
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_FLAGS "-fopenmp")

#Configuration du projet
project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})
#FIND_PACKAGE(BLAS)
#FIND_PACKAGE(LAPACK)
#set(MUMPS_DIR /opt/MUMPS_4.10.0)
#set(MPI_DIR /usr/lib/openmpi)
#set(BLAS /opt/OpenBLAS)
#set(LAPACK /opt/lapack)
#set(SLAPACK /opt/scalapack-2.0.1)
#set(BLACS /opt/BLACS/LIB)
#set(METIS /opt/ParMetis-3.1.1)
#set(SCOTCH /opt/scotch_5.1.11_esmumps)

#find_package (Boost COMPONENTS mpi REQUIRED)

link_directories(
  ../${CMAKE_BUILD_TYPE}
  /opt/openmpi-1.6.3-gnu/lib/
  #$ENV{HOME}/usr/lib/
  $ENV{HOME}/scratch/projects/sparselib_1_7/lib
  $ENV{HOME}/scratch/opt/mumps/intel/
  /opt/intel/mkl/lib/intel64/
  /opt/pt-scotch-5.1.12-gnu/lib/
  $ENV{HOME}/scratch/opt/boost_1_52_0/stage/lib
)

include_directories(../include $ENV{HOME}/scratch/opt/boost_1_52_0)

add_executable(
	test_abcd
	src/main.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
    test_abcd
    abcd 

    mkl_lapack95_lp64
    mkl_blas95_lp64
    mkl_blacs_openmpi_lp64
    mkl_scalapack_lp64
    mkl_blacs_openmpi_lp64
    mkl_intel_lp64
    mkl_gf_lp64 
    #mkl_gnu_thread
    mkl_sequential
    mkl_core 
    #pthread
    mkl_blacs_lp64
    mkl_sequential
    #mkl_gnu_thread
    mkl_core
    #pthread
    #pthread

  
    esmumps
    scotch
    scotcherr
    ptscotch
    
    dmumps
    mumps_common

    sparse
    spblas
    mv

    mpi
    mpi_f90
    mpi_f77
    mpi_cxx
    dl

    boost_mpi
    boost_serialization
)
