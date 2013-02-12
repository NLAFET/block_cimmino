cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_FC_COMPILER mpif90)
#SET(CMAKE_CXX_FLAGS "-std=c++11")

project(abcd CXX C Fortran)

include(abcdCmake.in)

find_package (Boost COMPONENTS mpi serialization REQUIRED)

# the library
#set(CMAKE_BUILD_TYPE "Debug")
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE}/lib)

#add_definitions(-DCOMPLEX="std::complex<double>")

include_directories(
  include
  ${Boost_INCLUDE_DIRS}
  ${SPARSELIB_DIR}/include
  ${MV_DIR}/include
  ${MUMPS_INCLUDE_DIR}
  )

file(
    GLOB_RECURSE
    source_files
    src/*.f
    src/*.cpp
)


add_library(abcd STATIC ${source_files})

add_subdirectory(test)

