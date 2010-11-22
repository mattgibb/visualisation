# Sets up CxxTest
INCLUDE(CxxTest.cmake)

SET(DATA_ROOT "${PROJECT_SOURCE_DIR}/test/data")
set(TEMP "${PROJECT_BINARY_DIR}/Testing/Temporary")

# Directory for test IO
FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")

# Hello world CxxTest
# CXXTEST_ADD_TEST(MyTestSuite_runner MyTestSuite_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MyTestSuite.h)

# itkVoxImageIO unit test
CXXTEST_ADD_TEST(itkVoxImageIOTest_runner itkVoxImageIOTest_runner.cpp )
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

SET(TEMP ${ITKTesting_BINARY_DIR}/Temporary)
SET(IO_TESTS ${CXX_TEST_PATH}/itkIOTests)

ADD_TEST(itkRawImageIOTest ${IO_TESTS} itkRawImageIOTest
                           ${TEMP}/raw1.raw ${TEMP}/raw2.raw
)

# TARGET_LINK_LIBRARIES(unittest_foo foo) # as needed

# enable dashboard scripting
# uses CTestConfig.cmake
# INCLUDE(CTest)
