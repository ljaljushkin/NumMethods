cmake_minimum_required(VERSION 2.8)
project(tbb)

include_directories("${CMAKE_SOURCE_DIR}/include")

set(SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/mmio.cpp
    ${CMAKE_SOURCE_DIR}/src/util.cpp)

set(HEADERS
    ${CMAKE_SOURCE_DIR}/include/timer.h
    ${CMAKE_SOURCE_DIR}/include/mmio.h
    ${CMAKE_SOURCE_DIR}/include/util.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

add_executable(solution
    ${SOURCES}
    ${HEADERS})
