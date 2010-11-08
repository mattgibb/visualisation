#ifndef __itkVoxImageIO_h
#define __itkVoxImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include "itkImageIOBase.h"
#include "itkImageRegion.h"
#include "itkPixelTraits.h"
#include "itkByteSwapper.h"
#include "itkVersion.h"
#include <string>

namespace itk
{
/** \class VoxImageIO
 *
 * \brief Read and write .vox ASCII images
 *
 * This class reads and writes 3D .vox images, compatible with Tarantula.
 *
 * \sa ImageFileReader
 *
 * \ingroup IOFilters
 */

class ITK_EXPORT VoxImageIO:public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef VoxImageIO           Self;
  typedef ImageIOBase          Superclass;
  typedef SmartPointer< Self > Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VoxImageIO, ImageIOBase);

  /** Pixel typedef support. Used to declare pixel type in filters
   * or other operations. */
  typedef char PixelType;

  /** this type is used in case the pixel has several components */
  typedef PixelTraits< PixelType >::ValueType ComponentType;

  /** Helper class to swap bytes when necessary */
  typedef ByteSwapper< ComponentType > ByteSwapperType;

  /** If the data is in the tail end of the file, you want to
   * explicitly set the header size. */
  void SetHeaderSize(unsigned long size);

  unsigned long GetHeaderSize();

  /** The different types of ImageIO's can support data of varying
   * dimensionality. For example, some file formats are strictly 2D
   * while others can support 2D, 3D, or even n-D. This method returns
   * true/false as to whether the ImageIO can support the dimension
   * indicated. */
  virtual bool SupportsDimension(unsigned long dim)
  { return ( dim == 3 ); }

  /*-------- This part of the interface deals with reading data. ------ */

  /** Determine the file type. Returns true if this ImageIOBase can read the
   * file specified. Always returns false because we don't want to use
   * this reader unless absolutely sure (i.e., manual ImageIO creation). */
  virtual bool CanReadFile(const char *);
  
  /** Read the first two lines for size and spacing */
  virtual void ReadImageInformation();
  
  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void *buffer);
  
  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Returns true if this ImageIO can write the specified file.
   * False is only returned when the file name is not specified. Otherwise
   * true is always returned. */
  virtual bool CanWriteFile(const char *);

  virtual void WriteImageInformation(void);

  /** Writes the data to disk from the memory buffer provided. */
  virtual void Write(const void *buffer);

protected:
  VoxImageIO();
  ~VoxImageIO();
  void PrintSelf(std::ostream & os, Indent indent) const;

  void OpenFileForReading(std::ifstream & is);

  void OpenFileForWriting(std::ofstream & os);

private:
  VoxImageIO(const Self &);     //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  std::string m_InternalFileName;

  bool           m_ManualHeaderSize;
  unsigned long  m_HeaderSize;
};
}
#endif
