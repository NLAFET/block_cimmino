cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpic++)
SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_FLAGS "-fopenmp")

#Configuration du projet
project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

find_package (Boost COMPONENTS mpi REQUIRED)

message(STATUS "ABCD_PROJECT_DIR=${ABCD_PROJECT_DIR}")
message(STATUS "ABCD_PROJECT_BINARY_DIR=${ABCD_PROJECT_BINARY_DIR}")

link_directories(
    ${ABCD_PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}
    ${ABCD_PROJECT_DIR}/lib
    ${Boost_LIBRARY_DIRS}
    ${LIB_DIRS}
)

include_directories(../include ${Boost_INCLUDE_DIRS})

add_executable(
        test_abcd
        src/main.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
    test_abcd
    abcd 
    patoh
    ${Boost_LIBRARIES}
    ${LIBS}
)
