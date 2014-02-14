CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( DEC )


#set(DEC_INCLUDE_INSTALL_DIR ${DEC_INCLUDE_INSTALL_DIR} CACHE PATH "The directory where we install the header files (optional)")
#
## set the internal install path for header files which depends on wether the user modifiable
## DEC_INCLUDE_INSTALL_DIR has been set by the user or not.
#if(DEC_INCLUDE_INSTALL_DIR)
#    set(INCLUDE_INSTALL_DIR
#        ${DEC_INCLUDE_INSTALL_DIR}
#        CACHE INTERNAL
#        "The directory where we install the header files (internal)"
#        )
#else()
#    set(INCLUDE_INSTALL_DIR
#        "${CMAKE_INSTALL_PREFIX}/include/DEC"
#        CACHE INTERNAL
#        "The directory where we install the header files (internal)"
#        )
#    message(INCLUDE_INSTALL_DIR)
#endif()








SET(CMAKE_RUNTIME_OUTPUT_DIR ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIR ${PROJECT_BINARY_DIR}/lib)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include $ENV{DEC_PATH})
INCLUDE_DIRECTORIES(/usr/local/include/eigen3)

ADD_DEFINITIONS(-std=c++11)
ADD_DEFINITIONS(-Wall)
ADD_DEFINITIONS(-Werror)

SET(DEC_SOURCES
    src/trianglemesh.cpp
    )
SET(DEC_HEADERS
    include/simplex.hpp
    include/simplicialComplex.hpp
    include/io.hpp
    include/types.hpp
    include/index.hpp
    include/dec.hpp
    include/dec_overloads.hpp
    include/util.hpp
    include/geometry.hpp
    include/render.hpp
    include/advection.hpp
    include/trianglemesh.h
    #include/kdtree.hpp
    )

ADD_CUSTOM_TARGET(headers SOURCES ${DEC_HEADERS})
ADD_LIBRARY(dec STATIC ${DEC_SOURCES} ${DEC_HEADERS})

#SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
#SET (CMAKE_CXX_COMPILER             "/usr/bin/g++")
INSTALL(FILES ${DEC_HEADERS} DESTINATION /usr/local/include/DEC)


ENABLE_TESTING()
ADD_SUBDIRECTORY(tests EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(tools)
ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)




##These lines are blatantly sotlen from eigen cmake file...
#string(TOLOWER "${CMAKE_GENERATOR}" cmake_generator_tolower)
#if(cmake_generator_tolower MATCHES "makefile")
#    message(STATUS "Some things you can do now:")
#    message(STATUS "--------------+--------------------------------------------------------------")
#    message(STATUS "Command       |   Description")
#    message(STATUS "--------------+--------------------------------------------------------------")
#    message(STATUS "make install  | Install to ${CMAKE_INSTALL_PREFIX}. To change that:")
#    message(STATUS "              |     cmake . -DCMAKE_INSTALL_PREFIX=yourpath")
#    message(STATUS "              |   DEC headers will then be installed to:")
#    message(STATUS "              |     ${INCLUDE_INSTALL_DIR}")
#    message(STATUS "              |   To install Eigen headers to a separate location, do:")
#    message(STATUS "              |     cmake . -DDEC_INCLUDE_INSTALL_DIR=yourpath")
#    message(STATUS "make doc      | Generate the API documentation, requires Doxygen & LaTeX")
#    message(STATUS "make check    | Build and run the unit-tests. Read this page:")
#    message(STATUS "              |   http://eigen.tuxfamily.org/index.php?title=Tests")
#    message(STATUS "make blas     | Build BLAS library (not the same thing as Eigen)")
#    message(STATUS "--------------+--------------------------------------------------------------")
#else()
#    message(STATUS "To build/run the unit tests, read this page:")
#    message(STATUS "  http://eigen.tuxfamily.org/index.php?title=Tests")
#endif()
#
#message(STATUS "")

