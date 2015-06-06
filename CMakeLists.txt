cmake_minimum_required(VERSION 3.2)
project(satellite_propagation_project)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")


include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE SOURCE_FILES "src/*.c")
file(GLOB TEST_FILES "tests/*.c" "tests/*.h")

add_library(satellite_propagation SHARED ${SOURCE_FILES})

add_executable(tests ${TEST_FILES} tests/main.c)
target_link_libraries(tests LINK_PUBLIC satellite_propagation)