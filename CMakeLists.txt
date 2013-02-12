cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_FC_COMPILER mpif90)
#SET(CMAKE_CXX_FLAGS "-std=c++11")

project(abcd CXX C Fortran)

#find_package (Boost COMPONENTS mpi serialization unit_test_framework REQUIRED)

# the library
#set(CMAKE_BUILD_TYPE "Debug")
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE}/lib)

#add_definitions(-DCOMPLEX="std::complex<double>")

include_directories(
  include
  #${Boost_INCLUDE_DIRS}
  #$ENV{HOME}/deploy/include/
  #$ENV{HOME}/usr/include/
  $ENV{HOME}/scratch/opt/mumps/intel/
  $ENV{HOME}/scratch/projects/sparselib_1_7/include/
  $ENV{HOME}/scratch/projects/sparselib_1_7/mv/include/
  $ENV{HOME}/scratch/opt/boost_1_52_0)

file(
    GLOB_RECURSE
    source_files
    src/*.f
    src/*.cpp
    include/*
)


add_library(abcd STATIC ${source_files})

add_subdirectory(test)
