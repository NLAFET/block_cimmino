cmake_minimum_required(VERSION 2.6)

project(Test_abcd)

if(DEBUG MATCHES ON)
  add_definitions(-DDEBUG)
endif()

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

if(NOT PATOH MATCHES OFF)
  set(LIBS ${LIBS} patoh)
endif()

set(LIBS ${LIBS}
  dmumps
  mumps_common

  ${Boost_LIBRARIES}
  ${BLAS_LAPACK_SCALAPACK_LIBS}
  ${MPI_LIBRARIES}

  sparse
  spblas
  mv
  )

link_directories(
  ${ABCD_PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}
  ${ABCD_PROJECT_DIR}/lib
  ${ABCD_PROJECT_DIR}/lib/sparselib/lib/
  ${ABCD_PROJECT_DIR}/lib/mumps/lib/
  ${Boost_LIBRARY_DIRS}
  ${MUMPS_LIBRARY_DIR}
  ${BLAS_LAPACK_SCALAPACK_DIRS}
  ${MPI_LIB_DIR}
  ${LIB_DIRS}
  )

include_directories(../include ${Boost_INCLUDE_DIRS})


add_executable(
  cpp_example
  ${ABCD_PROJECT_DIR}/example/example.cpp
  )


###Configuration de l'édition de liens
target_link_libraries(
  cpp_example
  abcd 
  ${LIBS}
  )

add_executable(
  abcd_run
  src/main.cpp
)


#Configuration de l'édition de liens
target_link_libraries(
  abcd_run
  abcd 
  ${LIBS}
)


################################
# Testing
################################
if(TEST MATCHES ON)

  add_executable(
          test_abcd_c_interface
          src/c_exemple.c
          src/c_interface.cc
  )


  #Configuration de l'édition de liens
  target_link_libraries(
    test_abcd_c_interface
    abcd 
    ${LIBS}
    )

  add_executable(
    basic_test
    src/basic.cpp
    )

  target_link_libraries(
    basic_test
    abcd 
    ${LIBS}
    )

  # Locate GTest

  set(GTEST_LIBRARY $ENV{HOME}/.local/lib/libgtest.a)
  set(GTEST_MAIN_LIBRARY $ENV{HOME}/.local/lib/libgtest_main.a)
  set(GTEST_INCLUDE_DIR $ENV{HOME}/.local/include/)

  set(GMOCK_LIBRARY $ENV{HOME}/.local/lib/libgmock.a)
  set(GMOCK_MAIN_LIBRARY $ENV{HOME}/.local/lib/libgmock_main.a)
  set(GMOCK_INCLUDE_DIR $ENV{HOME}/.local/include/)

  find_package(GTest REQUIRED)
  find_package(GMock REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})

  # Link runTests with what we want to test and the GTest and pthread library
  add_executable(runTests src/test.cpp)
  target_link_libraries(runTests
    ${GTEST_BOTH_LIBRARIES}
    ${GMOCK_BOTH_LIBRARIES}
    abcd 
    ${LIBS}
    )

  # Link runTests with what we want to test and the GTest and pthread library
  add_executable(runParallelTests src/test_para.cpp)
  target_link_libraries(runParallelTests
    ${GTEST_BOTH_LIBRARIES}
    ${GMOCK_BOTH_LIBRARIES}
    abcd 
    ${LIBS}
    )
  enable_testing()

endif()
