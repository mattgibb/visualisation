# Sets up CxxTest
INCLUDE(CxxTest.cmake)

SET(DATA_ROOT "${PROJECT_SOURCE_DIR}/test/data")
set(TEMP "${PROJECT_BINARY_DIR}/Testing/Temporary")

# Hello world CxxTest
CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)

# convert_image integration test
ADD_TEST(convert_image
         ${CONVERTERS_BINARY_DIR}/bin/convert_image
         ${DATA_ROOT}/input.vox
         ${TEMP}/output.vox
        )
        
ADD_TEST(convert_image_output
         diff
         ${DATA_ROOT}/input.vox
         ${TEMP}/output.vox
         )

ADD_EXECUTABLE(ImageReadWrite ImageReadWrite.cxx)
TARGET_LINK_LIBRARIES(ImageReadWrite ITKIO)

ADD_TEST(plugin_read_image
         ${CURRENT_BINARY_DIR}/ImageReadWrite
         ${DATA_ROOT}/input.vox
         ${TEMP}/read_write_image.mhd
        )

ADD_EXECUTABLE(ImageIO ImageIO.cxx)
TARGET_LINK_LIBRARIES(ImageIO ITKIO)

ADD_TEST(plugin_write_read_write
         ${CURRENT_BINARY_DIR}/ImageIO
         ${TEMP}/random_output.vox
         ${TEMP}/random_output2.mhd
        )
        
# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
