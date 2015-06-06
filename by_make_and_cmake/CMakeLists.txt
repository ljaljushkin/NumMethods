cmake_minimum_required(VERSION 2.8)
project(serial)

include_directories("${CMAKE_SOURCE_DIR}/include")

set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/mmio.cpp
    ${CMAKE_SOURCE_DIR}/src/util.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(solution ${SOURCE_FILES})