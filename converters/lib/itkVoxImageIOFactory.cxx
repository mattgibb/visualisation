#ifndef __itkVoxImageIOFactory_cxx
#define __itkVoxImageIOFactory_cxx

#include "itkVoxImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkVoxImageIO.h"
#include "itkVersion.h"

namespace itk
{
VoxImageIOFactory::VoxImageIOFactory()
{
  this->RegisterOverride( "itkImageIOBase",
                          "itkVoxImageIO",
                          "Vox Image IO",
                          1,
                          CreateObjectFunction< VoxImageIO >::New() );
}

VoxImageIOFactory::~VoxImageIOFactory()
{}

const char *
VoxImageIOFactory::GetITKSourceVersion() const
{
  return ITK_SOURCE_VERSION;
}

const char *
VoxImageIOFactory::GetDescription() const
{
  return "Vox ImageIO Factory, allows the reading and writing of Tarantula .vox images into insight";
}
} // end namespace itk
#endif
