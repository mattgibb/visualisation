#ifndef __itkVoxImageIOFactory_h
#define __itkVoxImageIOFactory_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
/** \class VoxImageIOFactory
 * \brief Create instances of VoxImageIO objects using an object factory.
 */
{
class ITK_EXPORT VoxImageIOFactory:public ObjectFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef VoxImageIOFactory Self;
  typedef ObjectFactoryBase                            Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion(void) const;

  virtual const char * GetDescription(void) const;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Not sure if we need this code */
  // static VoxImageIOFactory * FactoryNew() { return new VoxImageIOFactory; }
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(VoxImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    ObjectFactoryBase::RegisterFactory( Self::New() );
  }

protected:
  VoxImageIOFactory();
  ~VoxImageIOFactory();
private:
  VoxImageIOFactory(const Self &); //purposely not implemented
  void operator=(const Self &);    //purposely not implemented
};
} // namespace itk
#endif
