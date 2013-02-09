CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( DEC )


SET(CMAKE_RUNTIME_OUTPUT_DIR ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include $ENV{EIGEN_PATH})

ADD_DEFINITIONS(-std=c++11)
ADD_DEFINITIONS(-Wall)
#ADD_DEFINITIONS(-Werror)

SET(DEC_HEADERS
    include/simplex.hpp
    include/simplicialComplex.hpp
    include/io.hpp
    include/types.hpp
    include/index.hpp
    include/dec.hpp
    include/dec_overloads.hpp
    include/util.hpp
    include/render.hpp
    )

ADD_CUSTOM_TARGET(headers SOURCES ${DEC_HEADERS})

SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
#SET (CMAKE_CXX_COMPILER             "/usr/bin/g++")
INSTALL(FILES ${DEC_HEADERS} DESTINATION /usr/local/include/DEC)


ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(tools)

