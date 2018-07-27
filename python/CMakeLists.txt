cmake_minimum_required(VERSION 2.6)

##Build Type
set(BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

##Configuration file
include(../abcdCmake.in)

##Project
project(abcdpy)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

IF("${ABCD_PROJECT_DIR}" STREQUAL "")
    set(ABCD_PROJECT_DIR ${ABCD})
ENDIF()

##Environment file
IF(NOT ENV_INCLUDED)
    include(../env.in)
ENDIF()

#convert ${ABCD} to absolute
get_filename_component(ABCD_DIR "${ABCD}" ABSOLUTE)
add_subdirectory("${CMAKE_SOURCE_DIR}/python/numpy")

##Library abcdpy
add_library(abcdpy STATIC
  python_interface.cpp
)
set_target_properties(abcdpy PROPERTIES PREFIX "")
target_link_libraries(abcdpy
  abcd
  boost_numpy
  ${LIBS}
  )
