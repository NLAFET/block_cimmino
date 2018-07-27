cmake_minimum_required(VERSION 2.6)

##Build Type
set(BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

##Configuration file
include(../abcdCmake.in)

##Project
project(example CXX C Fortran)
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

IF ("${ABCD}" STREQUAL "")
    message(FATAL_ERROR "You need to provide the path to abcd root as: 'cmake .. -DABCD=/path/to/abcd'")
ENDIF()
IF("${ABCD_PROJECT_DIR}" STREQUAL "")
    set(ABCD_PROJECT_DIR ${ABCD})
ENDIF()
set(LIB_DIRS
    ${ABCD_PROJECT_DIR}/build/lib
    ${LIB_DIRS}
)

##Environment file
#IF(NOT ENV_INCLUDED)
    include(../env.in)
#ENDIF()

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
