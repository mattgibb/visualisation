#include "itkVoxPlugin.h"
#include "itkVoxImageIOFactory.h"

/**
 * Routine that is called when the shared library is loaded by
 * itk::ObjectFactoryBase::LoadDynamicFactories().
 *
 * itkLoad() is C (not C++) function.
 */
itk::ObjectFactoryBase* itkLoad()
{
  static itk::VoxImageIOFactory::Pointer f = itk::VoxImageIOFactory::New();
  return f;
}