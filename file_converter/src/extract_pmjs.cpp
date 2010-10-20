//
// Creates a set of pmjs for a free-running Purkinje pmj file and a surface definition file.
//
#include <sstream>

#include <boost/program_options.hpp>

#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkConfigure.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkSmartPointer.h"

using namespace boost;

namespace po = boost::program_options;

/**
 * Extracts candidate pmjs from the Purkinje system (terminal nodes)
 */
void get_pkj_pmjs(const vtkSmartPointer<vtkUnstructuredGrid>& pkj_grid,
                  std::vector<vtkIdType>& pmjs);

int main(int argc, char *argv[])
{
  //Process command line 
  //set up possible options
  po::options_description desc("Command line arguments");
  desc.add_options()
    ("help", "produce help message")
    ("surface", po::value<std::string>(), "Endocardial surface unstructured grid file")
    ("pkj", po::value<std::string>(), "Purkinje input file")
    ("output", po::value<std::string>(), "Output file name")
    ;

  //Parse command line
  po::variables_map vm; 
  po::store(po::parse_command_line(argc, argv, desc), vm); 
  po::notify(vm);

  if (vm.count("help")) 
  { 
    std::cout << desc << "\n"; 
    return EXIT_FAILURE;
  }

  //read pkj file
  vtkXMLUnstructuredGridReader* pkj_reader = vtkXMLUnstructuredGridReader::New();
  pkj_reader->SetFileName(vm["pkj"].as<std::string>().c_str());
  pkj_reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> pkj_grid = pkj_reader->GetOutput();
  vtkPoints* pkj_points = pkj_grid->GetPoints();

  std::vector<vtkIdType> pmjs;
  get_pkj_pmjs(pkj_grid, pmjs);

  //read surface file
  vtkXMLUnstructuredGridReader* surface_reader = vtkXMLUnstructuredGridReader::New();
  surface_reader->SetFileName(vm["surface"].as<std::string>().c_str());
  surface_reader->Update();
  vtkUnstructuredGrid* surface_grid = surface_reader->GetOutput();
  vtkPoints* surface_points = surface_grid->GetPoints();  
  
  //Setup point locator
  vtkPointLocator *point_locator = vtkPointLocator::New();
  point_locator->SetDataSet(surface_grid);
  point_locator->BuildLocator();

  std::ofstream pmjs_file(vm["output"].as<std::string>().c_str());
  if (!pmjs_file)
  {
    std::cerr << "Failed to open .pmj file" << std::endl;
    abort();
  }

  //Create Vertices
  for (std::vector<vtkIdType>::iterator iter = pmjs.begin();
       iter != pmjs.end();
       ++iter)
  {
    vtkIdType pkj_pt_id = *iter;
 
    pmjs_file << pkj_pt_id << " ";

    //Find closest surface point
    vtkIdType surface_pt_id = point_locator->FindClosestPoint(pkj_points->GetPoint(pkj_pt_id));

    /* std::cout << "point: " << pkj_points->GetPoint(pkj_pt_id)[0] << " "
              << pkj_points->GetPoint(pkj_pt_id)[1] << " "
              << pkj_points->GetPoint(pkj_pt_id)[2] << " surface: "
              << surface_points->GetPoint(surface_pt_id)[0] << " "
              << surface_points->GetPoint(surface_pt_id)[1] << " "
              << surface_points->GetPoint(surface_pt_id)[2] << std::endl;*/
    
    pmjs_file << surface_pt_id << std::endl;
  }
  
  pmjs_file.close();

  //clean up
  pkj_reader->Delete();
  surface_reader->Delete();
  point_locator->Delete();
}

//
//
//
void get_pkj_pmjs(const vtkSmartPointer<vtkUnstructuredGrid>& pkj_grid,
                  std::vector<vtkIdType>& pmjs)
{
  //Loop over points, initializing to zeros
  vtkSmartPointer<vtkPoints> pkj_points = pkj_grid->GetPoints();

  std::vector<unsigned> candidates(pkj_points->GetNumberOfPoints());
  
  //loop over VTK lines, increment each index
  vtkSmartPointer<vtkIdList> points = vtkSmartPointer<vtkIdList>::New();
  points->SetNumberOfIds(2);

  for (vtkIdType i = 0; i < pkj_grid->GetNumberOfCells(); ++i)
  {
    pkj_grid->GetCellPoints(i, points);

    ++candidates[points->GetId(0)];
    ++candidates[points->GetId(1)];
  }

  for (vtkIdType i = 0; i < candidates.size(); ++i)
  {
    if(candidates[i] == 1)
    {
      pmjs.push_back(i);
    }
  }
}
