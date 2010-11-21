# Sets up CxxTest
INCLUDE(CxxTest.cmake)

# Hello world test
# CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)

SET(DATA_ROOT "${PROJECT_SOURCE_DIR}/test/data")
set(TEMP "${PROJECT_BINARY_DIR}/Testing/Temporary")

# convert_image integration test
ADD_TEST(convert_image
         ${EXECUTABLE_OUTPUT_PATH}/convert_image
         ${DATA_ROOT}/input.vox
         ${TEMP}/output.vox
        )

ADD_EXECUTABLE(ImageReadWrite ImageReadWrite.cxx)
TARGET_LINK_LIBRARIES(ImageReadWrite ITKIO)

ADD_TEST(read_image
         ${EXECUTABLE_OUTPUT_PATH}/ImageReadWrite
         ${DATA_ROOT}/input.vox
         ${TEMP}/read_write_image.mhd
        )

# TARGET_LINK_LIBRARIES(unittest_foo foo) # as needed

# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
