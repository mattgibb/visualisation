cmake_minimum_required(VERSION 2.8)

PROJECT(FILE_CONVERTER)

# Find Boost and add include dirs
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

# Find Boost and add include dirs
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${USE_VTK_FILE})

# add project root to include path
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# MESSAGE( STATUS "HOME: $ENV{HOME}" )

  
# cheeky little shortcut if a group of things are needed.
# sets the varable ITKLibs to the rest of the arguments.
# SET(ITKLibs ITKCommon ITKBasicFilters ITKIO ITKNumerics)
SET(VTKLibs vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon)

# INCLUDE( CMakePrintSystemInformation )

ADD_EXECUTABLE(tetgen_to_vtk tetgen_to_vtk.cpp )
TARGET_LINK_LIBRARIES(tetgen_to_vtk ${Boost_LIBRARIES} vtkCommon vtkIO)

# GET_TARGET_PROPERTY(LINK_FLAGS REGISTER LINK_FLAGS)
