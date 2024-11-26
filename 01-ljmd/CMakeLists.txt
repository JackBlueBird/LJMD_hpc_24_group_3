cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(LJMD VERSION 1.0)

# for MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# Use configure_file to set LJMD version
configure_file(LJMDConfig.h.in LJMDConfig.h @ONLY)

# add the header files to the list of source files
set(HEADER_FILES src/mdsys.h src/read_input.h src/utilities.h src/ekin.h src/force.h src/verlet.h src/output.h)
# Add separated files into a library
add_library(mdlib SHARED src/mdsys.c src/read_input.c src/utilities.c src/ekin.c src/force.c src/verlet.c src/output.c ${HEADER_FILES})
#target_compile_options(mdlib PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer )
target_compile_options(mdlib PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer -D_MPI)
target_link_libraries(mdlib m ${MPI_LIBRARIES})
# Include source file directory
target_include_directories(mdlib PUBLIC ./src)

# add the executable
add_executable(ljmd-serial.x src/main.c)

# for MPI
#target_link_libraries(ljmd-serial.x ${MPI_LIBRARIES})

target_compile_options(ljmd-serial.x PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
# link math.h and mdlib
target_link_libraries(ljmd-serial.x m mdlib ${MPI_LIBRARIES})
# Use target_include_directories to include ${PROJECT_BINARY_DIR}
target_include_directories(ljmd-serial.x PUBLIC "${PROJECT_BINARY_DIR}" ./src)


# Testing
option(ENABLE_TESTING "Enable building unit tests" ON)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  #----------------------------------------
  add_executable(test_verlet gtest/tests/test_verlet.cpp)
  target_compile_options(test_verlet PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/src)
  add_test(NAME Verlet COMMAND test_verlet)

  set_target_properties(test_verlet
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )

  #----------------------------------------
  add_executable(test_utils gtest/tests/test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Utils COMMAND test_utils)

  set_target_properties(test_utils
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )

  #----------------------------------------
  add_executable(test_read_inp gtest/tests/test_read_inp.cpp)
  target_include_directories(test_read_inp PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_read_inp PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME ReadInputData
    COMMAND sh -c "${CMAKE_CURRENT_SOURCE_DIR}/build/gtest-bin/test_read_inp < argon_3.inp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/gtest/examples"
  )

  set_target_properties(test_read_inp
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )

  #----------------------------------------
  add_executable(test_read_rest gtest/tests/test_read_rest.cpp)
  target_include_directories(test_read_rest PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_read_rest PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME ReadRestartData
    COMMAND sh -c "${CMAKE_CURRENT_SOURCE_DIR}/build/gtest-bin/test_read_rest < argon_3.inp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/gtest/examples"
  )

  set_target_properties(test_read_rest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )

  #----------------------------------------
  add_executable(test_ekin gtest/tests/test_ekin.cpp)
  target_include_directories(test_ekin PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Ekin COMMAND test_ekin)

  set_target_properties(test_ekin
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )

  #----------------------------------------
  add_executable(test_force gtest/tests/test_force.cpp)
  target_include_directories(test_force PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Force COMMAND test_force)

  set_target_properties(test_force
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gtest-bin"
  )


endif()
