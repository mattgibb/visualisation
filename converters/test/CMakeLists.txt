# Sets up CxxTest
INCLUDE(CxxTest.cmake)

SET(DATA_ROOT "${PROJECT_SOURCE_DIR}/test/data")
set(TEMP "${PROJECT_BINARY_DIR}/Testing/Temporary")

# Hello world CxxTest
CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)

ADD_EXECUTABLE(ImageReadWrite ImageReadWrite.cxx)
TARGET_LINK_LIBRARIES(ImageReadWrite ITKIO)

# convert_image integration test
ADD_TEST(convert_image
         ${CONVERTERS_BINARY_DIR}/bin/convert_image
         ${DATA_ROOT}/input.vox
         ${TEMP}/output.vox
        )

ADD_TEST(plugin_read_image
         ${EXECUTABLE_OUTPUT_PATH}/ImageReadWrite
         ${DATA_ROOT}/input.vox
         ${TEMP}/read_write_image.mhd
        )

# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
