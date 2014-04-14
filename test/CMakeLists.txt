cmake_minimum_required(VERSION 2.6)

project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

if(NOT PATOH MATCHES OFF)
    set(LIBS ${LIBS} patoh)
endif()

set(LIBS ${LIBS}
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES}
    ${BLAS_LAPACK_SCALAPACK_LIBS}
    dmumps
    mumps_common

    sparse
    spblas
    mv
)

link_directories(
    ${ABCD_PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}
    ${ABCD_PROJECT_DIR}/lib
    ${ABCD_PROJECT_DIR}/lib/sparselib/lib/
    ${Boost_LIBRARY_DIRS}
    ${MUMPS_LIBRARY_DIR}
    ${BLAS_LAPACK_SCALAPACK_DIRS}
    ${MPI_LIB_DIR}
    ${LIB_DIRS}
)

include_directories(../include ${Boost_INCLUDE_DIRS})


add_executable(
        cpp_example
        src/example.cpp
)


###Configuration de l'édition de liens
target_link_libraries(
    cpp_example
    abcd 
    ${LIBS}
)


################################
# Testing
################################
if(test MATCHES ON)
    add_executable(
            test_abcd
            src/main.cpp
    )


    #Configuration de l'édition de liens
    target_link_libraries(
        test_abcd
        abcd 
        ${LIBS}
    )

    add_executable(
            test_abcd_c_interface
            src/c_test.c
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

    set(GTEST_LIBRARY $ENV{HOME}/.local/lib/libgtest.a)
    set(GTEST_MAIN_LIBRARY $ENV{HOME}/.local/lib/libgtest_main.a)
    set(GTEST_INCLUDE_DIR $ENV{HOME}/.local/include/)

    set(GMOCK_LIBRARY $ENV{HOME}/.local/lib/libgmock.a)
    set(GMOCK_MAIN_LIBRARY $ENV{HOME}/.local/lib/libgmock_main.a)
    set(GMOCK_INCLUDE_DIR $ENV{HOME}/.local/include/)

  # Locate GTest
  find_package(GTest REQUIRED)
  find_package(GMock REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})

  # Link runTests with what we want to test and the GTest and pthread library
  add_executable(runTests src/test.cpp)
  message("${GMOCK_BOTH_LIBRARIES}")
  target_link_libraries(runTests
      ${GTEST_BOTH_LIBRARIES}
      ${GMOCK_BOTH_LIBRARIES}
      abcd 
      ${LIBS}
  )

  enable_testing()

  ##############
  # Unit Tests
  ##############
  #add_executable(runUnitTests src/test_utils.cpp)

  ## Standard linking to gtest stuff.
  #target_link_libraries(runUnitTests
      #${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES}
      #abcd 
      #${Boost_LIBRARIES}
      #${LIBS}
  #)

  # Extra linking for the project.
  #target_link_libraries(runUnitTests project1_lib)

  # This is so you can do 'make test' to see all your tests run, instead of
  # manually running the executable runUnitTests to see those specific tests.
  #add_test(NAME test COMMAND runTests)
endif()
