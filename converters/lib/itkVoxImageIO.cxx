#ifndef __itkVoxImageIO_cxx
#define __itkVoxImageIO_cxx
#include "itkVoxImageIO.h"

namespace itk
{
VoxImageIO::VoxImageIO():
  ImageIOBase()
{
  this->SetNumberOfComponents(1);
  this->SetPixelTypeInfo( typeid( char ) );
  this->SetNumberOfDimensions(3);

  for ( unsigned int idx = 0; idx < 3; ++idx )
    {
    m_Spacing.insert(m_Spacing.begin() + idx, 1.0);
    m_Origin.insert(m_Origin.begin() + idx, 0.0);
    }
  
  m_HeaderSize = 0;
  m_ManualHeaderSize = false;
  
  // Left over from short reader
  m_ByteOrder = ImageIOBase::BigEndian;
  m_FileType = ASCII;
  
  // set origin as in itkPNGImageIO
  m_Origin[0] = 0.0;
  m_Origin[1] = 0.0;
  m_Origin[2] = 0.0;
}

VoxImageIO::~VoxImageIO()
{}

void VoxImageIO::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "I'm a vox reader!" << std::endl;
}

unsigned long VoxImageIO::GetHeaderSize()
{
  std::ifstream file;

  if ( m_FileName == "" )
    {
    itkExceptionMacro(<< "A FileName must be specified.");
    }
    
  // not sure about order of initialisation here, so commented it out.
  // if ( !m_ManualHeaderSize )
  //   {
  //   this->ComputeStrides();
  // 
  // make sure we figure out a filename to open
  this->OpenFileForReading(file);
  
  //   // Get the size of the header from the size of the image
  //   file.seekg(0, std::ios::end);
  // 
  //   m_HeaderSize = (unsigned long)( (unsigned long)file.tellg()
  //                                   - (unsigned long)m_Strides[4] * 2 );
  //   }
    
  file.seekg(0, std::ios::beg);
      
  char dummy[300];
  file.getline(dummy,300).getline(dummy,300);
  
  if ( file.tellg() == -1 )
    {
    itkExceptionMacro(<< "The first couple of vox file lines are longer than 300 characters.")
    }
  
  m_HeaderSize = file.tellg();
  this->Modified();
    // if ((unsigned long)file.tellg() != m_HeaderSize)
    //   {
    //   itkExceptionMacro(<< "The vox file is not formatted properly.");
    //   }
  file.close();
  
  return m_HeaderSize;
}

void VoxImageIO::OpenFileForReading(std::ifstream & is)
{
  if ( m_FileName == "" )
    {
    itkExceptionMacro(<< "A FileName must be specified.");
    }

  // Close file from any previous image
  if ( is.is_open() )
    {
    is.close();
    }

  // Open the new file
  itkDebugMacro(<< "Initialize: opening file " << m_FileName);
  is.open(m_FileName.c_str(), std::ios::in);
  if ( is.fail() )
    {
    itkExceptionMacro(<< "Could not open file: " << m_FileName);
    }
}

void VoxImageIO::OpenFileForWriting(std::ofstream & os)
{
  if ( m_FileName == "" )
    {
    itkExceptionMacro(<< "A FileName must be specified.");
    }

  // Close file from any previous image
  if ( os.is_open() )
    {
    os.close();
    }

  // Open the new file
  itkDebugMacro(<< "Initialize: opening file " << m_FileName);
  os.open(m_FileName.c_str(), std::ios::out);
  if ( os.fail() )
    {
    itkExceptionMacro(<< "Could not open file: " << m_FileName);
    }
}

void VoxImageIO::SetHeaderSize(unsigned long size)
{
  if ( size != m_HeaderSize )
    {
    m_HeaderSize = size;
    this->Modified();
    }
  m_ManualHeaderSize = true;
}

bool VoxImageIO::CanReadFile(const char *name)
{
  std::string filename = name;
  
  if (filename == "" )
    {
      return false;
    }
    
  std::string::size_type voxPos = filename.rfind(".vox");
  if ( ( voxPos != std::string::npos )
       && ( voxPos == filename.length() - 4 ) )
    {
    return true;
    }
  
  return false;
}

void VoxImageIO::ReadImageInformation()
{
  std::ifstream file;
  
  // Make sure m_Dimensions and m_Spacing are empty
  m_Dimensions.clear();
  m_Spacing.clear();
  
  // Open the file
  this->OpenFileForReading(file);
  
  // initialise m_Dimensions and m_Spacing
  SizeValueType size;
  std::cout << "Before m_Dimensions, file.tellg(): " << file.tellg() << std::endl;
  for (unsigned int i=0; i<3; i++)
  {
    file >> size;
    std::cout << "In m_Dimensions, i = " << i << ", file.tellg(): " << file.tellg() << std::endl;
    m_Dimensions.push_back(size);
  }
  
  double spacing;
  for (unsigned int i=0; i<3; i++)
  {
    file >> spacing;
    std::cout << "In m_Spacing, i = " << i << ", file.tellg(): " << file.tellg() << std::endl;
    m_Spacing.push_back(spacing);
  }
  
  for (unsigned int i=0; i<3; i++)
  {
    std::cout << "m_Dimensions[" << i << "]: " << m_Dimensions[i] << std::endl;
    std::cout << "m_Spacing[" << i << "]: " << m_Spacing[i] << std::endl;
  }
  
  std::cout << "m_Dimensions.size(): " << m_Dimensions.size() << std::endl;
  std::cout << "m_Spacing.size(): " << m_Spacing.size() << std::endl;
  
  file.close();
}

void VoxImageIO::Read(void *buffer)
{
  std::ifstream file;

  // Open the file
  this->OpenFileForReading(file);
  this->ComputeStrides();
  
  // Offset into file
  unsigned long streamStart = this->GetHeaderSize();
  file.seekg( (long)streamStart, std::ios::beg );
  if ( file.fail() )
    {
    itkExceptionMacro(<< "File seek failed");
    }

  const unsigned long numberOfBytesToBeRead =
    static_cast< unsigned long >( this->GetImageSizeInBytes() );

  itkDebugMacro(<< "Reading " << numberOfBytesToBeRead << " bytes");

  this->ReadBufferAsASCII( file, buffer, this->GetComponentType(),
                           this->GetImageSizeInComponents() );

  itkDebugMacro(<< "Reading Done");

#define itkReadRawBytesAfterSwappingMacro(StrongType, WeakType)   \
  ( this->GetComponentType() == WeakType )                        \
    {                                                             \
    typedef ByteSwapper< StrongType > InternalByteSwapperType;    \
    if ( m_ByteOrder == LittleEndian )                            \
      {                                                           \
      InternalByteSwapperType::SwapRangeFromSystemToLittleEndian( \
        (StrongType *)buffer, this->GetImageSizeInComponents() ); \
      }                                                           \
    else if ( m_ByteOrder == BigEndian )                          \
      {                                                           \
      InternalByteSwapperType::SwapRangeFromSystemToBigEndian(    \
        (StrongType *)buffer, this->GetImageSizeInComponents() ); \
      }                                                           \
    }

  // Swap bytes if necessary
  if itkReadRawBytesAfterSwappingMacro(unsigned short, USHORT)
  else if itkReadRawBytesAfterSwappingMacro(short, SHORT)
  else if itkReadRawBytesAfterSwappingMacro(char, CHAR)
  else if itkReadRawBytesAfterSwappingMacro(unsigned char, UCHAR)
  else if itkReadRawBytesAfterSwappingMacro(unsigned int, UINT)
  else if itkReadRawBytesAfterSwappingMacro(int, INT)
  else if itkReadRawBytesAfterSwappingMacro(long, LONG)
  else if itkReadRawBytesAfterSwappingMacro(unsigned long, ULONG)
  else if itkReadRawBytesAfterSwappingMacro(float, FLOAT)
  else if itkReadRawBytesAfterSwappingMacro(double, DOUBLE)
}

bool VoxImageIO::CanWriteFile(const char *fname)
{
  std::string filename(fname);

  if ( filename == "" )
    {
    return false;
    }
  
  std::string::size_type voxPos = filename.rfind(".vox");
  if ( ( voxPos != std::string::npos )
       && ( voxPos == filename.length() - 4 ) )
    {
    return true;
    }
  
  
  return false;
}

void VoxImageIO::WriteImageInformation(void)
{}

void VoxImageIO::Write(const void *buffer)
{
  std::ofstream file;
  
  // Open the file
  this->OpenFileForWriting(file);
  
  // Set up for reading
  this->ComputeStrides();

  for(unsigned int i=0; i<m_Strides.size(); i++)
  {
    std::cout << "m_Strides[" << i << "]: " << m_Strides[i] << std::endl;
  }
  
  // Write the header
  file << this->GetDimensions(0) << " "
       << this->GetDimensions(1) << " "
       << this->GetDimensions(2) << "\n";
	file << this->GetSpacing(0)    << " "
	     << this->GetSpacing(1)    << " "
	     << this->GetSpacing(2)    << "\n";
  
  // Write the body
  this->WriteBufferAsVoxFormat( file, buffer, this->GetComponentType(), this->GetImageSizeInComponents() );
  
  std::cout << "this->GetImageSizeInBytes(): " << this->GetImageSizeInBytes() << std::endl;
  std::cout << "this->GetImageSizeInComponents(): " << this->GetImageSizeInComponents() << std::endl;
  std::cout << "this->GetImageSizeInPixels(): " << this->GetImageSizeInPixels() << std::endl;
  std::cout << "this->GetComponentType(): " << this->GetComponentType() << std::endl;
  std::cout << "this->GetComponentTypeAsString(this->GetComponentType()): " << this->GetComponentTypeAsString(this->GetComponentType()) << std::endl;
  
  file.close();
}

// similar to itk::ImageIOBase::WriteBuffer, but doesn't write 6 numbers per line.
template< class TComponent >
void WriteBufferOnePerLine(std::ostream & os, const TComponent *buffer, ImageIOBase::SizeType num)
{
  const TComponent *ptr = buffer;

  typedef typename itk::NumericTraits< TComponent >::PrintType PrintType;
  for ( ImageIOBase::SizeType i = 0; i < num; i++ )
    {
    os << PrintType(*ptr++) << "\n";
    }
}

// similar to itk::ImageIOBase::WriteBufferAsASCII, but uses WriteBufferOnePerLine instead of WriteBuffer.
// can't override ImageIOBase::WriteBufferAsASCII, as it's not virtual.
// Loads of boilerplate, if anyone has a better way of doing this, please do!
void VoxImageIO::WriteBufferAsVoxFormat(std::ostream & os, const void *buffer,
                                         IOComponentType ctype,
                                         ImageIOBase::SizeType numComp)
{
  switch ( ctype )
    {
    case UCHAR:
      {
      typedef const unsigned char *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;
    case CHAR:
      {
      typedef const char *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case USHORT:
      {
      typedef const unsigned short *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case SHORT:
      {
      typedef const short *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case UINT:
      {
      typedef const unsigned int *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case INT:
      {
      typedef const int *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case ULONG:
      {
      typedef const unsigned long *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case LONG:
      {
      typedef const long *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case FLOAT:
      {
      typedef const float *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    case DOUBLE:
      {
      typedef const double *Type;
      Type buf = reinterpret_cast< Type >( buffer );
      WriteBufferOnePerLine(os, buf, numComp);
      }
      break;

    default:
      break;
    }
}

} // namespace itk
#endif
