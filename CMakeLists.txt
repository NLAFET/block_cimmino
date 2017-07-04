cmake_minimum_required(VERSION 2.6)

set(BUILD_TYPE "Release")

if(DEBUG MATCHES ON)
    set(BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "-O0 --coverage")
    add_definitions(-DDBG)
endif()

include(abcdCmake.in)
project(abcd CXX C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package (Boost COMPONENTS mpi serialization REQUIRED)

# the library
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib)
set(ABCD_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(ABCD_PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR})

SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -m64 -z muldefs -std=c++11")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -m64 -z muldefs -std=c++11")
SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -m64 -z muldefs -std=c++11")

#add_definitions(-DCOMPLEX="std::complex<double>")

if(NOT PATOH MATCHES OFF)
    add_definitions(-DPATOH)
endif()

if(MUMPS_ES MATCHES OFF)
    add_definitions(-DNO_MUMPS_ES)
endif()

if(WIP MATCHES ON)
  add_definitions(-DWIP)
endif()

include_directories(
  ${ABCD_PROJECT_DIR}/include
  ${ABCD_PROJECT_DIR}/lib/sparselib/include
  ${ABCD_PROJECT_DIR}/lib/sparselib/mv/include
#  ${ABCD_PROJECT_DIR}/lib/mumps/include
  ${Boost_INCLUDE_DIRS}
  ${BOOST_ROOT}
  ${MUMPS_INCLUDE_DIR}
  ${MPI_INC_DIR}
  ${INC_DIRS}
  )

file(
    GLOB_RECURSE
    source_files
    src/*.f
    src/*.cpp
    src/*.c
)

add_library(abcd SHARED ${source_files})

add_subdirectory(test)

add_subdirectory(src/utils)
