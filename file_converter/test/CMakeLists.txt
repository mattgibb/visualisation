# Sets up CxxTest
INCLUDE(CxxTest.cmake)

CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)
# TARGET_LINK_LIBRARIES(unittest_foo foo) # as needed

# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
