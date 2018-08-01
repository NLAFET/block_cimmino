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

##Environment file
IF("${BOOST_ROOT}" STREQUAL "")
    message("BOOST_ROOT not specified. Compilation of Boost libraries mpi/serialization/python with ABCD sources.")
    set(BOOST_ROOT ${ABCD_PROJECT_DIR}/lib/boost/)
    set(Boost_INCLUDE_DIRS
        ${BOOST_ROOT}
        ${BOOST_ROOT}/boost
        ${BOOST_ROOT}/libs
    )
    set(Boost_LIBRARY_DIRS ${BOOST_ROOT}/../../build/lib)
    # check if libraries where already installed
    find_library(BOOSTMPI_FOUND NAMES libboost_mpi.a  libboost_serialization.a PATHS ${BOOST_LIBRARY_DIRS})
    find_library(BOOSTSER_FOUND NAMES libboost_serialization.a PATHS ${BOOST_LIBRARY_DIRS})
    find_library(BOOSTPYT_FOUND NAMES libboost_python.a PATHS ${BOOST_LIBRARY_DIRS})
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
    IF(NOT "${PYTHON_DEF}" STREQUAL "")
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

##Libs and Tests
# libabcd
file(
    GLOB_RECURSE
    source_files
    src/*.f
    src/*.cpp
    src/*.c
)
add_library(abcd STATIC ${source_files})
target_link_libraries(
    abcd
    ${LIBS}
)

# abcd_run
add_subdirectory(test)

# Python
IF(NOT "${PYTHON_DEF}" STREQUAL "")
    find_package (PythonLibs REQUIRED)
    if(PYTHONLIBS_FOUND)
        add_subdirectory(python)
    else()
        message("PythonLibs didn't work")
    endif()
endif()
