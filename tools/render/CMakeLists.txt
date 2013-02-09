
INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/include
    ./include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    )
SET(QTRENDER_SOURCES
    src/qtmain.cpp
    src/qtmainwindow.cpp
    src/glutil.cpp
    src/qtglwidget.cpp
    )
SET(QTRENDER_HEADERS_TOMOC
    include/qtmainwindow.h
    include/qtglwidget.h
    )
SET(QTRENDER_HEADERS
    ${QTRENDER_HEADERS_TOMOC}
    include/glutil.h
    )
SET(QTRENDER_RESOURCE
    shaders/shaders.qrc
    )
SET(QTRENDER_SHADERS
    shaders/shader.v.glsl
    shaders/shader.130.v.glsl
    shaders/shader.f.glsl
    shaders/shader.130.f.glsl
    shaders/shader.g.glsl
    shaders/faceshader.g.glsl
    shaders/edgeshader.g.glsl
    shaders/vertshader.g.glsl
    )
#This target has an obvious purpose :)
ADD_CUSTOM_TARGET(MAKE_SHADERS_VISIBLE_IN_QTCREATOR
    SOURCES ${QTRENDER_SHADERS})

FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)

QT4_ADD_RESOURCES(QTRENDER_RCC ${QTRENDER_RESOURCE})
MESSAGE(STATUS ${QTRENDER_RCC})
QT4_WRAP_CPP(QTRENDER_HEADERS_MOC
    ${QTRENDER_HEADERS_TOMOC})
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

ADD_EXECUTABLE(qtRender ${QTRENDER_SOURCES} ${QTRENDER_HEADERS} ${QTRENDER_RCC} ${QTRENDER_HEADERS_MOC} ${QTRENDER_SHADERS})
TARGET_LINK_LIBRARIES(qtRender ${QT_LIBRARIES} GLEW GL)
