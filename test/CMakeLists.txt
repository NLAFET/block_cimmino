cmake_minimum_required(VERSION 2.6)

##Configuration file
include(../abcdCmake.in)

##Project
IF("${PROJECT_NAME}" STREQUAL "")
    project(Test_abcd CXX)

    IF ("${ABCD}" STREQUAL "")
        message(FATAL_ERROR "You need to provide the path to abcd root as: 'cmake .. -DABCD=/path/to/abcd'")
    ENDIF()
    set(ABCD_PROJECT_DIR ${ABCD})

    ##Environment file
    IF(NOT ENV_INCLUDED)
        include(../env.in)
    ENDIF()
ENDIF()

# Executable
set(EXECUTABLE_OUTPUT_PATH ../)

##Executable example for C++
# Compiled main
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
if(NOT "${GTEST}" STREQUAL "")
#if(TEST MATCHES ON)
  # Locate GTest/GMock
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
			#    abcd
    ${LIBS}
    )

  # Link runTests with what we want to test and the GTest and pthread library
  add_executable(runParallelTests src/test_para.cpp)
  target_link_libraries(runParallelTests
    ${GTEST_BOTH_LIBRARIES}
    ${GMOCK_BOTH_LIBRARIES}
			#    abcd
    ${LIBS}
    )
  enable_testing()
endif()
