
INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/include
    ./include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    /home/mtao/hg/eigen/unsupported
    )



ADD_EXECUTABLE(voronoi src/main.cpp)
TARGET_LINK_LIBRARIES(voronoi dec)


