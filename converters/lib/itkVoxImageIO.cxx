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
  m_Origin[3] = 0.0;
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

  if ( !m_ManualHeaderSize )
    {
    this->ComputeStrides();

    // make sure we figure out a filename to open
    this->OpenFileForReading(file);

    // Get the size of the header from the size of the image
    file.seekg(0, std::ios::end);

    m_HeaderSize = (unsigned long)( (unsigned long)file.tellg()
                                    - (unsigned long)m_Strides[4] * 2 );
    }
    
    file.seekg(0, std::ios::beg);
    
    std::cout << "tellg() after seekg(0, std::ios::beg): " << file.tellg() << std::endl;
    
    char dummy[300];
    file.getline(dummy,300).getline(dummy,300);
    
    std::cout << "file.tellg(): " << file.tellg() << std::endl;
    std::cout << "m_HeaderSize: " << m_HeaderSize << std::endl;
    
    if ((unsigned long)file.tellg() != m_HeaderSize)
      {
      itkExceptionMacro(<< "The vox file is not formatted properly.");
      }
    
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
  
  // Open the file
  this->OpenFileForReading(file);
  
  // initialise m_Dimensions and m_Spacing
  SizeValueType size;
  for (unsigned int i=0; i<3; i++)
  {
    file >> size;
    m_Dimensions.push_back(size);
  }
  
  double spacing;
  for (unsigned int i=0; i<3; i++)
  {
    file >> spacing;
    m_Spacing.push_back(spacing);
  }
  
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

  // Write the header
  file << this->GetDimensions(0) << " "
       << this->GetDimensions(1) << " "
       << this->GetDimensions(2) << "\n";
	file << this->GetSpacing(0)    << " "
	     << this->GetSpacing(1)    << " "
	     << this->GetSpacing(2)    << "\n";

  // Write the body
  const char *buf = reinterpret_cast< const char * >( buffer );
  long num = this->GetImageSizeInBytes();
  
  for ( ImageIOBase::SizeType i = 0; i < num; i++ )
  {
    file << buf++ << "\n";
  }
  
  file.close();
}
} // namespace itk
#endif
