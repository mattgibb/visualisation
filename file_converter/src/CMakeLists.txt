# add current directory to include path
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# cheeky little shortcut if a group of things are needed.
# sets the varable VTKLibs to the rest of the arguments.
SET(VTKLibs vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon)

# tetgen_to_vtk target and libraries to link to
ADD_EXECUTABLE(tetgen_to_vtk tetgen_to_vtk.cpp )
TARGET_LINK_LIBRARIES(tetgen_to_vtk ${Boost_LIBRARIES} vtkCommon vtkIO)
