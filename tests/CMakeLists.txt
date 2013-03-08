ENABLE_TESTING()

INCLUDE_DIRECTORIES(. ../include)
ADD_EXECUTABLE(compileTest compile.cpp)
ADD_EXECUTABLE(ioTest io.cpp)
ADD_EXECUTABLE(operatorTest operators.cpp)

ADD_EXECUTABLE(d0p_3Test d0p_3.cpp)


ADD_TEST(compileTest ${CMAKE_CURRENT_BINARY_DIR}/compileTest)
ADD_TEST(ioTest ${CMAKE_CURRENT_BINARY_DIR}/ioTest)
ADD_TEST(operatorTest ${CMAKE_CURRENT_BINARY_DIR}/operatorTest)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS compileTest ioTest operatorTest)
