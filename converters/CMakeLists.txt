cmake_minimum_required(VERSION 2.8.2)

PROJECT(CONVERTERS)

# Find Boost and add include dirs
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/lib)

INCLUDE(CTest)

# Find libraries and add include dirs
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${USE_VTK_FILE})
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

if(NOT EXECUTABLE_OUTPUT_PATH)
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
endif(NOT EXECUTABLE_OUTPUT_PATH)

ADD_SUBDIRECTORY(lib)

# Execute ${PROJECT_SOURCE_DIR}/src/CMakeLists.txt, and put outputs in ${PROJECT_BUILD_DIR}/bin
ADD_SUBDIRECTORY(src bin)

# Testing stuff
OPTION(BUILD_TESTING "Set this to build and run the testing framework")
if(BUILD_TESTING)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(test)
endif(BUILD_TESTING)
