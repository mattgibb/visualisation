// Uses ITK to convert any image format ITK understands to any format ITK understands.
// If necessary e.g. for custom/proprietory file formats, can write an ImageIO subclass
// and register it with the IO factory.

#include <boost/program_options.hpp>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  // PARSE COMMAND LINE OPTIONS
  po::positional_options_description pd;
  pd.add("input", 1).add("output", 1);
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "Display this help message.")
    ("input", po::value<std::string>(), "Name of input file.  Alternative syntax to first positional parameter.")
    ("output", po::value<std::string>(), "name of output file. Alternative syntax to second positional parameter.")
  ;
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(desc).positional(pd).run(), vm);
  po::notify(vm);
 
  // if --help is specified, or the usage is wrong
  if ( vm.count("help") || !vm.count("input") || !vm.count("output") ) {
    cout << "\nUsage: " << argv[0] << " input output\n\n";
    cout << desc << "\n";
    return EXIT_FAILURE;
  }
  
  
  // ITK STUFF
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];

  reader->SetFileName( inputFilename  );
  writer->SetFileName( outputFilename );

  writer->SetInput( reader->GetOutput() );

  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
    
  return EXIT_SUCCESS;
  
  
}
