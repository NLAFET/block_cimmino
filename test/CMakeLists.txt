cmake_minimum_required(VERSION 2.6)

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)

#Configuration du projet
project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

find_package (Boost COMPONENTS mpi REQUIRED)

link_directories(../${CMAKE_BUILD_TYPE} ${Boost_LIBRARY_DIRS})
include_directories(../include ${Boost_INCLUDE_DIRS})

add_executable(
	test_abcd
	src/main.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
	test_abcd
	${Boost_LIBRARIES}
	abcd 
)