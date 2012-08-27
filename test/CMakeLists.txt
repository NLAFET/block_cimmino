cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

#Configuration du projet
project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

find_package (Boost COMPONENTS mpi unit_test_framework REQUIRED)


# Dependencies
set(PETSC_DIR /home/knuthy/clones/petsc-3.3-p2/)
set(BLAS_DIR /opt/OpenBLAS)

set(
  PETSC_INCLUDE_DIR 
  ${PETSC_DIR}/include 
  ${PETSC_DIR}/linux/include
  )
  
set(PETSC_LIBRARY_DIR ${PETSC_DIR}/linux/lib)

link_directories ( ../include ${Boost_LIBRARY_DIRS} )
include_directories ( ../${CMAKE_BUILD_TYPE}/lib/ ${Boost_INCLUDE_DIRS} )

# add_executable(
# 	test_abcd
# 	src/main.cpp
# )

add_executable(
  boost_test_abcd
  src/test.cpp
)
enable_testing()
add_test(boost_test_abcd boost_test_abcd)

#Configuration de l'édition de liens
# target_link_libraries(
# 	test_abcd
# 	abcd
# )
target_link_libraries(
  boost_test_abcd
  abcd ${Boost_LIBRARIES}
  mpi
#   blacs
#   lapack
#   blas
#   X11
#   dl
#   scalapack
#   parmetis
#   metis
#   petsc
#   dmumps
#   mumps_common
#   scalapack
#   blacs
#   pord
  )