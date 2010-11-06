# add current directory to include path
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# cheeky little shortcut if a group of things are needed.
# sets the varable VTKLibs to the rest of the arguments.
SET(VTKLibs vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon)

# targets and libraries to link to
ADD_EXECUTABLE(convert_carp_to_vtk convert_carp_to_vtk.cpp )
TARGET_LINK_LIBRARIES(convert_carp_to_vtk ${Boost_LIBRARIES} vtkCommon vtkIO)

ADD_EXECUTABLE(convert_chaste_to_vtk convert_chaste_to_vtk.cpp )
TARGET_LINK_LIBRARIES(convert_chaste_to_vtk ${Boost_LIBRARIES} vtkCommon vtkIO)
