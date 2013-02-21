
INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/include
    ./include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    /home/mtao/hg/eigen/unsupported
    )
SET(LIBDECRENDER_SOURCES
    src/qtmainwindow.cpp
    src/glutil.cpp
    src/qtglwidget.cpp
    )
SET(LIBDECRENDER_HEADERS_TOMOC
    include/qtmainwindow.h
    include/qtglwidget.h
    )
SET(LIBDECRENDER_HEADERS
    ${LIBDECRENDER_HEADERS_TOMOC}
    include/glutil.h
    include/packages.h
    )
SET(LIBDECRENDER_RESOURCE
    shaders/shaders.qrc
    )
SET(LIBDECRENDER_SHADERS
    shaders/shader.v.glsl
    shaders/shader.130.v.glsl
    shaders/shader.f.glsl
    shaders/shader.130.f.glsl
    shaders/noneshader.130.v.glsl
    shaders/shader.g.glsl
    shaders/faceshader.g.glsl
    shaders/edgeshader.g.glsl
    shaders/vertshader.g.glsl
    )
#This target has an obvious purpose :)
ADD_CUSTOM_TARGET(MAKE_SHADERS_VISIBLE_IN_QTCREATOR
    SOURCES ${LIBDECRENDER_SHADERS})

FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)

QT4_ADD_RESOURCES(LIBDECRENDER_RCC ${LIBDECRENDER_RESOURCE})
MESSAGE(STATUS ${LIBDECRENDER_RCC})
QT4_WRAP_CPP(LIBDECRENDER_HEADERS_MOC
    ${LIBDECRENDER_HEADERS_TOMOC})
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

ADD_LIBRARY(decrender STATIC ${LIBDECRENDER_SOURCES} ${LIBDECRENDER_HEADERS} ${LIBDECRENDER_RCC} ${LIBDECRENDER_HEADERS_MOC} ${LIBDECRENDER_SHADERS})

TARGET_LINK_LIBRARIES(decrender ${QT_LIBRARIES} GLEW GL)


ADD_EXECUTABLE(qtRender src/qtmain.cpp)
TARGET_LINK_LIBRARIES(qtRender ${QT_LIBRARIES} decrender)


ADD_EXECUTABLE(formExample src/formexample.cpp)
TARGET_LINK_LIBRARIES(formExample ${QT_LIBRARIES} decrender)
#TARGET_LINK_LIBRARIES(qtRender ${QT_LIBRARIES} GLEW GL arpack)
