/*
 * Converts Chaste output vtk files to paraview timeseries compatible pvd output
 */
#include <cmath>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

using namespace boost;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  //Process command line 
  //set up possible options
  po::options_description desc("Command line arguments");
  desc.add_options()
    ("help", "produce help message")
    ("input", po::value<std::string>(), "name of input file")
    ("output_root", po::value<std::string>(), "The root of the output file, will have XXX.vtu appended")
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

  if (!vm.count("input")) 
  { 
    std::cout << "Input was not set.\n";
    return EXIT_FAILURE;
  } 

  if (!vm.count("output_root")) 
  { 
    std::cout << "Output was not set.\n";
    return EXIT_FAILURE;
  } 

  //Load Chaste file.
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(vm["input"].as<std::string>().c_str());
  reader->Update();

  vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();

  vtkSmartPointer<vtkPointData> point_data = grid->GetPointData();

  vtkSmartPointer<vtkUnstructuredGrid> output_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  output_grid->SetPoints(grid->GetPoints());
  output_grid->SetCells(grid->GetCellTypesArray(),
                        grid->GetCellLocationsArray(),
                        grid->GetCells());

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  //We need someway to flag which nodes are intracellular. This data doesn't seem to be outputted by 
  //Chaste, so we retain the intial v_0 data to threshold on. Ideally this should be a deep copy to 
  //avoid naming issues.
  vtkSmartPointer<vtkDataArray> v0_array = point_data->GetArray(0);
  v0_array->SetName("V0");
  output_grid->GetPointData()->AddArray(v0_array);

  //Loop over arrays, there should be two arrays for each time level
  // in the order Vm, phi_e. Note we start at timestep=1 to avoid deep copy issues.
  for (int i = 2; i < point_data->GetNumberOfArrays()/2; ++i)
  {
    std::cout << "Converting " << i << " of " << point_data->GetNumberOfArrays()/2 << std::endl;

    vtkSmartPointer<vtkDataArray> v_array = point_data->GetArray(2*i);
    vtkSmartPointer<vtkDataArray> phi_array = point_data->GetArray(2*i + 1);

    v_array->SetName("Vm");
    phi_array->SetName("phi_e");

    output_grid->GetPointData()->AddArray(v_array);
    output_grid->GetPointData()->AddArray(phi_array);

    writer->SetInput(output_grid);
    
    //Setup the file name for output. Format needs to be output_rootXXX.vtu
    std::stringstream num;
    num.width(3);
    num <<  std::setfill('0') << std::right << i;
   
    std::cout << num.str() << std::endl;

    writer->SetFileName((vm["output_root"].as<std::string>() + num.str() + ".vtu").c_str());
    writer->Write();
  }
}
