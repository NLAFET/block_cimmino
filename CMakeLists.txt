cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_FC_COMPILER mpif90)

project(abcd CXX C Fortran)

find_package (Boost COMPONENTS mpi serialization unit_test_framework REQUIRED)

# Dependencies
set(EIGEN_DIR /usr/local/include/eigen3)

# the library
set(CMAKE_BUILD_TYPE "Debug")
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE}/lib)


include_directories(include ${Boost_INCLUDE_DIRS} ${EIGEN_DIR})

file(
	GLOB_RECURSE
	source_files
	src/*.f
	src/*
	include/*
)

add_library(abcd STATIC ${source_files})

add_subdirectory(test)
