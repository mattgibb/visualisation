# Sets up CxxTest
INCLUDE(CxxTest.cmake)

# Hello world test
CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)

# convert_image integration test
ADD_TEST(convert_image
         ${PROJECT_BINARY_DIR}/src/convert_image
         ${CMAKE_CURRENT_SOURCE_DIR}/data/input.vox
         ${CMAKE_CURRENT_SOURCE_DIR}/data/output.vox
        )

# TARGET_LINK_LIBRARIES(unittest_foo foo) # as needed

# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
