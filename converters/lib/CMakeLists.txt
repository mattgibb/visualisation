SET(VOXIO_SRCS itkVoxImageIO.cxx itkVoxImageIOFactory.cxx)

ADD_LIBRARY(VOXIO ${VOXIO_SRCS})
TARGET_LINK_LIBRARIES(VOXIO ITKIO)

# Shared library that when placed in ITK_AUTOLOAD_PATH, will add VoxImageIO
# as an ImageIOFactory without recompiling for any ITK application.
ADD_LIBRARY(VoxPlugin SHARED itkVoxPlugin.cxx)
TARGET_LINK_LIBRARIES(VoxPlugin VOXIO)
