
cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

add_library(mdlib src/utils.c src/verlet.c src/read_input.c src/ekin.c src/force.c)

option(ENABLE_TESTING "Enable building unit tests" OFF)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  add_executable(test_verlet tests/test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils tests/test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Utils COMMAND test_utils)

  # add_executable(test_read tests/test_read.cpp)
  # target_include_directories(test_read PRIVATE ${CMAKE_SOURCE_DIR}/src)
  # target_link_libraries(test_read PRIVATE GTest::GTestMain GTest::GTest mdlib)
  # add_test(NAME Readinput
  #   COMMAND sh -c "${CMAKE_CURRENT_SOURCE_DIR}/build/test_read < argon_3.inp"
  #   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples"
  # )

  add_executable(test_read_inp tests/test_read_inp.cpp)
  target_include_directories(test_read_inp PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_read_inp PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME ReadInputData
    COMMAND sh -c "${CMAKE_CURRENT_SOURCE_DIR}/build/test_read_inp < argon_3.inp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples"
  )

  add_executable(test_read_rest tests/test_read_rest.cpp)
  target_include_directories(test_read_rest PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_read_rest PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME ReadRestartData
    COMMAND sh -c "${CMAKE_CURRENT_SOURCE_DIR}/build/test_read_rest < argon_3.inp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples"
  )

  add_executable(test_ekin tests/test_ekin.cpp)
  target_include_directories(test_ekin PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Ekin COMMAND test_ekin)

  add_executable(test_force tests/test_force.cpp)
  target_include_directories(test_force PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Force COMMAND test_force)

endif()
