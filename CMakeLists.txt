cmake_minimum_required(VERSION 3.2)
project(satellite_propagation)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES satellite_propagation/main.c)
add_executable(satellite_propagation ${SOURCE_FILES} satellite_propagation/rotation_matrix.h satellite_propagation/rotation_matrix.c satellite_propagation/constants.h satellite_propagation/constants.c satellite_propagation/date_converters/date_converters.h satellite_propagation/date_converters/date_converters.c satellite_propagation/precession_matrix.h satellite_propagation/precession_matrix.c satellite_propagation/matrix_operations.h satellite_propagation/matrix_operations.c)