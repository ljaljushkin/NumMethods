cmake_minimum_required(VERSION 2.8)

#SET(CMAKE_C_COMPILER /opt/intel/composer_xe_2015.3.187/bin/intel64/icс)
#SET(CMAKE_CXX_COMPILER /opt/intel/composer_xe_2015.3.187/bin/intel64/icpc)
#
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

project(openmp)

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
