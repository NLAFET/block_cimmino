cmake_minimum_required(VERSION 2.6)

##Build Type
set(BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

##Configuration file
include(../abcdCmake.in)

IF(NOT "${PROJECT_NAME}" STREQUAL "")
    message(FATAL_ERROR "example should be compiled alone: 'cd example/; mkdir build; cd build; cmake .. -DABCD=/path/to/abcd; make'")
ENDIF()
IF ("${ABCD}" STREQUAL "")
    message(FATAL_ERROR "You need to provide the path to abcd root as: 'cmake .. -DABCD=/path/to/abcd'")
ENDIF()

##Project
project(example CXX C Fortran)
set(ABCD_PROJECT_DIR ${ABCD})
set(LIB_DIRS
    ${ABCD_PROJECT_DIR}/build/lib
    ${LIB_DIRS}
)

include(../env.in)

##Executable example for C++
# Compiled main
add_executable(
  example
  example.cpp
  )
# Configuration de l'édition de liens
target_link_libraries(
  example
  abcd
  ${LIBS}
  )

##Executable example for C
# Compiled main
add_executable(
  example_c
  example.c
  )
# Configuration de l'édition de liens
target_link_libraries(
  example_c
  abcd
  ${LIBS}
  stdc++
  m
  )
