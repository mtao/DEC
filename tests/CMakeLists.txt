ENABLE_TESTING()

INCLUDE_DIRECTORIES(. ../include)
ADD_EXECUTABLE(compileTest compile.cpp)
ADD_TEST(compileTest ${CMAKE_CURRENT_BINARY_DIR}/compileTest)

