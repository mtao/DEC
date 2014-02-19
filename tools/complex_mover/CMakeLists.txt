cmake_minimum_required(VERSION 2.8.9)
PROJECT(complex_mover)

SET(CMAKE_RUNTIME_OUTPUT_DIR ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)

ADD_DEFINITIONS(-std=c++11)
ADD_DEFINITIONS(-Werror)
ADD_DEFINITIONS(-Wall)

# Include dirs
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include /opt/intel/tbb/include/ /home/mtao/git/DEC/include)


LINK_DIRECTORIES(/usr/local/lib64 /opt/intel/tbb/lib/intel64/gcc4.4)

# Source files
SET(complex_mover_SOURCES
    src/main.cpp
    )

SET(complex_mover_HEADERS
    )

# Executable to build
ADD_EXECUTABLE(complex_mover
    ${complex_mover_SOURCES}
    ${complex_mover_HEADERS}
    )

TARGET_LINK_LIBRARIES(complex_mover openvdb tbb Half dec)
