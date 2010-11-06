cmake_minimum_required(VERSION 2.8.2)

PROJECT(FILE_CONVERTER)

# Find Boost and add include dirs
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

# Find libraries and add include dirs
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${USE_VTK_FILE})
FIND_PACKAGE(ITK REQUIRED)

ADD_SUBDIRECTORY(src)


OPTION(BUILD_TESTING "Set this to build and run the testing framework")
if(BUILD_TESTING)
  ADD_SUBDIRECTORY(test)
endif(BUILD_TESTING)
