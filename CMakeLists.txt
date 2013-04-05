cmake_minimum_required(VERSION 2.6)


project(abcd CXX C Fortran)

include(abcdCmake.in)

find_package (Boost COMPONENTS mpi serialization REQUIRED)

# the library
#set(CMAKE_BUILD_TYPE "Debug")
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE}/lib)
set(ABCD_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(ABCD_PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR})

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

