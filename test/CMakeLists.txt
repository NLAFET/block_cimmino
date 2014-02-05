cmake_minimum_required(VERSION 2.6)

project(Test_abcd)

set(EXECUTABLE_OUTPUT_PATH ../${CMAKE_BUILD_TYPE})

#find_package (Boost COMPONENTS mpi REQUIRED)

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

add_executable(
        basic_test
        src/basic.cpp
)

#Configuration de l'édition de liens
target_link_libraries(
    test_abcd
    abcd 
    ${Boost_LIBRARIES}
    ${LIBS}
)


################################
# Testing
################################
if(test MATCHES ON)
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
      ${Boost_LIBRARIES}
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
