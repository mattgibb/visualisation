#ifndef __itkVoxPlugin_h
#define __itkVoxPlugin_h

#include "itkObjectFactoryBase.h"

#ifdef WIN32
#ifdef VoxPlugin_EXPORTS
#define VoxPlugin_EXPORT __declspec(dllexport)
#else
#define VoxPlugin_EXPORT __declspec(dllimport)
#endif
#else
#define VoxPlugin_EXPORT 
#endif

/**
 * Routine that is called when the shared library is loaded by
 * itk::ObjectFactoryBase::LoadDynamicFactories().
 *
 * itkLoad() is C (not C++) function.
 */
extern "C" {
VoxPlugin_EXPORT itk::ObjectFactoryBase* itkLoad();
} 
#endif
