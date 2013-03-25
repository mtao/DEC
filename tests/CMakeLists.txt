ENABLE_TESTING()

INCLUDE_DIRECTORIES(. ../include)
ADD_EXECUTABLE(compileTest compile.cpp)
ADD_EXECUTABLE(ioTest io.cpp)
ADD_EXECUTABLE(operatorTest operators.cpp)


ADD_TEST(compileTest ${CMAKE_CURRENT_BINARY_DIR}/compileTest)
ADD_TEST(ioTest ${CMAKE_CURRENT_BINARY_DIR}/ioTest)
ADD_TEST(operatorTest ${CMAKE_CURRENT_BINARY_DIR}/operatorTest)

SET(DECOPS d0p d1p d0d d1d h0p h1p h2p h0d h1d h2d)

FOREACH(DECOP ${DECOPS})
    SET(TNAME ${DECOP} "Test")
    MESSAGE("[" ${DECOP} "] [" ${TNAME} "] [" ${DECOP}_3.cpp "]")
    #ADD_EXECUTABLE(${TNAME} ${DECOP}_3.cpp)
    #ADD_TEST(${TNAME} ${CMAKE_CURRENT_BINARY_DIR}/${TNAME})
ENDFOREACH(DECOP)






add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS compileTest ioTest operatorTest)
