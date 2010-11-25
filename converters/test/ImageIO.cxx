#include <fstream>
#include "itkRandomImageSource.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

#include "cxxtest/TestSuite.h"

using namespace std;

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc != 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile " << std::endl;
    return EXIT_FAILURE;
  }

  typedef unsigned short PixelType;
  typedef itk::Image< PixelType,3 > ImageType;
  typedef itk::ImageRegionConstIterator< ImageType > ImageIteratorType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  // set filenames
  string vox1(argv[1]);
  string vox2(argv[2]);
  
  // Create a source object (in this case a random image generator).
  // The source object is templated on the output type.
  //
  ImageType::SizeValueType size[3];

  size[0]=32; size[1]=16; size[2]=8;

  itk::RandomImageSource<ImageType>::Pointer random;
  random = itk::RandomImageSource<ImageType>::New();
  random->SetMin(0);
  random->SetMax(24680);
  random->SetSize(size);

  // Create a mapper (in this case a writer). A mapper
  // is templated on the input type.
  //

  // Write out the image
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(random->GetOutput());
  writer->SetFileName(vox1.c_str());

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while writing the image " << vox1 << std::endl;
    std::cerr << excp << std::endl;
    }

  // Create a source object (in this case a reader)
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(vox1.c_str());

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while reading the image " << vox1 << std::endl;
    std::cerr << excp << std::endl;
    }


  // Compare pixel by pixel in memory
  ImageIteratorType it( reader->GetOutput(),
                        reader->GetOutput()->GetBufferedRegion() );

  ImageIteratorType ot( random->GetOutput(),
                        random->GetOutput()->GetBufferedRegion() );

  it.GoToBegin();
  ot.GoToBegin();
  while( !it.IsAtEnd() )
    {
    const PixelType iv = it.Get();
    const PixelType ov = ot.Get();
    if( iv != ov )
      {
      cerr << "Error in read/write of pixel " << it.GetIndex() << endl;
      cerr << "Read value  is : " << iv << endl;
      cerr << "it should be   : " << ov << endl;
      cerr << "Test FAILED ! " << endl;
      return EXIT_FAILURE;
      }
    ++it;
    ++ot;
    }
  
  writer->SetInput(reader->GetOutput());
  writer->SetFileName(vox2.c_str());
  writer->SetInput(reader->GetOutput());

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error while writing the image " << argv[2] << std::endl;
    std::cerr << excp << std::endl;
    }


  std::cerr << "Test PASSED ! " << std::endl;
  return EXIT_SUCCESS;
};
