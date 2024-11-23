cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(LJMD VERSION 1.0)

# Use configure_file to set LJMD version
configure_file(LJMDConfig.h.in LJMDConfig.h @ONLY)

# add the header files to the list of source files
set(HEADER_FILES src/mdsys.h src/read_input.h src/utilities.h src/ekin.h src/force.h src/velverlet.h src/output.h)
# Add separated files into a library
add_library(mdlib SHARED src/mdsys.c src/read_input.c src/utilities.c src/ekin.c src/force.c src/velverlet.c src/output.c ${HEADER_FILES})
# Include source file directory
target_include_directories(mdlib PUBLIC ./src)

# add the executable
add_executable(ljmd-serial.x src/main.c)
# link math.h and mdlib
target_link_libraries(ljmd-serial.x m mdlib)
# Use target_include_directories to include ${PROJECT_BINARY_DIR}
target_include_directories(ljmd-serial.x PUBLIC "${PROJECT_BINARY_DIR}" ./src)

