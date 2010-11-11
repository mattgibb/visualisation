cmake_minimum_required(VERSION 2.8.2)

PROJECT(FILE_CONVERTER)

# Find Boost and add include dirs
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/lib)

# Find libraries and add include dirs
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${USE_VTK_FILE})
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_SUBDIRECTORY(lib)

ADD_SUBDIRECTORY(src)

# Testing stuff
OPTION(BUILD_TESTING "Set this to build and run the testing framework")
if(BUILD_TESTING)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(test)
endif(BUILD_TESTING)
