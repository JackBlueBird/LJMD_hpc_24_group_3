cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(LJMD VERSION 1.0)

# for MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Use configure_file to set LJMD version
configure_file(LJMDConfig.h.in LJMDConfig.h @ONLY)

# add the header files to the list of source files
set(HEADER_FILES src/mdsys.h src/read_input.h src/utilities.h src/ekin.h src/force.h src/verlet.h src/output.h)
# Add separated files into a library
add_library(mdlib SHARED src/mdsys.c src/read_input.c src/utilities.c src/ekin.c src/force.c src/verlet.c src/output.c ${HEADER_FILES})
target_compile_options(mdlib PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer )
#target_compile_options(mdlib PUBLIC -Wall -g -O3 -ffast-math -fomit-frame-pointer -D_MPI)
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
