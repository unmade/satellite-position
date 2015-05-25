cmake_minimum_required(VERSION 3.2)
project(satellite_propagation)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES satellite_propagation/main.c)
add_executable(satellite_propagation ${SOURCE_FILES} satellite_propagation/rotation_matrix.h satellite_propagation/rotation_matrix.c satellite_propagation/constants.h satellite_propagation/constants.c satellite_propagation/date_converters/date_converters.h satellite_propagation/date_converters/date_converters.c satellite_propagation/precession.h satellite_propagation/precession.c satellite_propagation/matrix_operations.h satellite_propagation/matrix_operations.c satellite_propagation/nutation.h satellite_propagation/nutation.c satellite_propagation/coordinates_converters.h satellite_propagation/coordinates_converters.c satellite_propagation/forces/gravitational_potential.h satellite_propagation/forces/gravitational_potential.c)