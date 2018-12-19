cmake_minimum_required(VERSION 2.6)

##Build Type
set(BUILD_TYPE "Release")
if(DEBUG MATCHES ON)
    set(BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "-O0 --coverage")
    add_definitions(-DDBG)
    add_definitions(-DDEBUG)
endif()
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

##Configuration file
include(abcdCmake.in)

##Project
project(abcd CXX C Fortran)
set(ABCD_PROJECT_DIR ${PROJECT_SOURCE_DIR})
set(ABCD_PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR})

IF(DEFINED ENV{PYTHON_INC})
  set(PYTHON_INC $ENV{PYTHON_INC})
  message(STATUS "Add python include path : ${PYTHON_INC}")
  include_directories(${PYTHON_INC})
  set(PYTHON_DEF 1)
ELSE()
  message(FATAL_ERROR " You have to set up PYTHON_INC env variable with\
  the path to python include")
endif()

message(STATUS "PYTHON_DEF var = ${PYTHON_DEF}")

##Environment file
IF("${BOOST_ROOT}" STREQUAL "")
    message("BOOST_ROOT not specified. Compilation of Boost libraries mpi/serialization/python with ABCD sources.")
    set(BOOST_ROOT ${ABCD_PROJECT_DIR}/lib/boost/)
    message(STATUS "ABCD project dir = ${ABCD_PROJECT_DIR}")
    set(BOOST_INCLUDE_DIRS
        ${BOOST_ROOT}
        ${BOOST_ROOT}/boost
        ${BOOST_ROOT}/libs
    )
   #set(BOOST_LIBRARY_DIRS ${BOOST_ROOT}/../../build/lib)
    set(BOOST_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/lib)
    message(STATUS "BOOST_LIBRARY_DIRS = ${BOOST_LIBRARY_DIRS}")

################################################################################
#                  check if libraries were already installed                   #
################################################################################
    find_library(BOOSTMPI_FOUND NAMES libboost_mpi.a  libboost_serialization.a PATHS ${BOOST_LIBRARY_DIRS} NO_DEFAULT_PATH)
    message(STATUS "BOOSTMPI_FOUND = ${BOOSTMPI_FOUND}")
    IF("${BOOSTMPI_FOUND}" STREQUAL "BOOSTMPI_FOUND-NOTFOUND")
        # libboost_mpi
        file(
            GLOB_RECURSE
            src_boost
            lib/boost/libs/mpi/*.cpp
        )
        add_library(boost_mpi STATIC ${src_boost})
    ELSE()
        message("libboost_mpi already installed.")
    ENDIF()

    find_library(BOOSTSER_FOUND NAMES libboost_serialization.a PATHS ${BOOST_LIBRARY_DIRS} NO_DEFAULT_PATH)
    message(STATUS "BOOSTSER_FOUND = ${BOOSTSER_FOUND}")
    IF("${BOOSTSER_FOUND}" STREQUAL "BOOSTSER_FOUND-NOTFOUND")
        # libboost_serialization
        file(
          GLOB_RECURSE
          src_boost
          lib/boost/libs/serialization/*.cpp
          )
        add_library(boost_serialization STATIC ${src_boost})
    ELSE()
        message("libboost_serialization already installed.")
    ENDIF()

    IF(${PYTHON_DEF} EQUAL 1)
      find_library(BOOSTPYT_FOUND NAMES libboost_python.a PATHS ${BOOST_LIBRARY_DIRS} NO_DEFAULT_PATH)
      message(STATUS "BOOSTPYT_FOUND = ${BOOSTPYT_FOUND}")
      IF("${BOOSTPYT_FOUND}" STREQUAL "BOOSTPYT_FOUND-NOTFOUND")
        # libboost_python
        file(
          GLOB_RECURSE
          src_boost
          lib/boost/libs/python/*.cpp
          )
        add_library(boost_python STATIC ${src_boost})
      ELSE()
        message("libboost_python already installed.")
      ENDIF()
    ENDIF()
ELSE()
    message("BOOST_ROOT specified at ${BOOST_DIR}.")
    message("Checking Boost installation.")
    find_package (Boost COMPONENTS mpi serialization REQUIRED)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_NO_SYSTEM_PATHS ON)
ENDIF()

IF(NOT ENV_INCLUDED)
    include(env.in)
ENDIF()

#IF(DEFINED ENV{SPLDLT_LIB})
#  set(SPLDLT_LIB $ENV{SPLDLT_LIB} CACHE PATH "Location of the SpLDLT library")
#
#  find_library(SPLDLT_LIBRARY_FOUND NAMES spldlt PATHS ${SPLDLT_LIB})
# #message(STATUS "State of SPLDLT_LIBRARY_FOUND : ${SPLDLT_LIBRARY_FOUND}")
#
#  IF(SPLDLT_LIBRARY_FOUND)
#    set(SPLDLT_LIBRARY ${SPLDLT_LIBRARY_FOUND})
#    set(LIBS ${LIBS} ${SPLDLT_LIBRARY})
#    message(STATUS "Set SPLDLT_LIBRARY : ${SPLDLT_LIBRARY_FOUND}")
#  ELSE()
#    message(WARNING "libspldlt.a cannot be found in ${SPLDLT_LIB}")
#  ENDIF()
#ENDIF()
#
#
#IF(DEFINED ENV{SPLDLT_INC})
#  set(SPLDLT_INC $ENV{SPLDLT_INC} CACHE PATH "Location of the SpLDLT include")
#  #TODO add check in the include path
#  include_directories(${SPLDLT_INC})
#ELSE()
#  set(WARNING "No SPLDLT_INC include folder path variable set in the env")
#ENDIF()
#
#IF(DEFINED ENV{SPRAL_LIB})
#  set(SPRAL_LIB $ENV{SPRAL_LIB} CACHE PATH "Location of the Spral library")
#
#  find_library(SPRAL_LIBRARY_FOUND NAMES spral PATHS ${SPRAL_LIB})
# #message(STATUS "State of SPRAL_LIBRARY_FOUND : ${SPRAL_LIBRARY_FOUND}")
#
#  IF(SPRAL_LIBRARY_FOUND)
#    set(SPRAL_LIBRARY ${SPRAL_LIBRARY_FOUND})
#    set(LIBS ${LIBS} ${SPRAL_LIBRARY})
#    message(STATUS "Set SPRAL_LIBRARY : ${SPRAL_LIBRARY_FOUND}")
#  ELSE()
#    message(WARNING "libspral.a cannot be found in ${SPRAL_LIB}")
#  ENDIF()
#ENDIF()
#
#IF(DEFINED ENV{SPRAL_INC})
#  set(SPRAL_INC $ENV{SPRAL_INC} CACHE PATH "Location of the Spral include")
#  #TODO add check in the include path
#  include_directories(${SPRAL_INC})
#ELSE()
#  set(WARNING "No SPRAL_INC include folder path variable set in the env")
#ENDIF()


##Libs and Tests
# libabcd
file(
    GLOB_RECURSE
    source_files
    src/*.f
    src/*.cpp
    src/*.c
)
#message(STATUS "Make lib with ${source_files}")
add_library(abcd STATIC ${source_files})
target_link_libraries(
    abcd
    ${LIBS}
)

# abcd_run
add_subdirectory(test)

# Python
#IF(NOT "${PYTHON_DEF}" STREQUAL "")
#    find_package (PythonLibs REQUIRED)
#    if(PYTHONLIBS_FOUND)
#        add_subdirectory(python)
#    else()
#      message(ERROR "PythonLibs didn't work")
#    endif()
#endif()
