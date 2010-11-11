/*
 * converts meshalyzer mesh files to vtk unstructured data
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
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConfigure.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

using namespace boost;

namespace po = boost::program_options;

/**
 * Reads a .pts file and adds the points to the given grid
 */
void create_points(const std::string& file_name, 
                   vtkSmartPointer<vtkUnstructuredGrid> grid);

/**
 * Reads a .tris file and adds the triangles to the given grid
 */
void create_tris(const std::string& file_name,
                 vtkSmartPointer<vtkUnstructuredGrid> grid);

/**
 * Reads a .elem file and adds the triangles to the given grid
 */
void create_elem(const std::string& file_name,
                 vtkSmartPointer<vtkUnstructuredGrid> grid);

/**
 * Reads a .cnnx file and adds the lines to the given grid;
 */
void create_lines(const std::string& file_name,
                  vtkSmartPointer<vtkUnstructuredGrid> grid);

/**
 * Reads a .lon file and adds the vectors to the given elements;
 */

void create_vectors(const std::string& file_name,
                  vtkSmartPointer<vtkUnstructuredGrid> grid);


int main(int argc, char *argv[])
{
  //Process command line 
  //set up possible options
  po::options_description desc("Command line arguments");
  desc.add_options()
    ("help", "produce help message")
    ("pts", po::value<std::string>(), "name of input pts file")
    ("tris", po::value<std::string>(), "name of input tris file")
    ("elem", po::value<std::string>(), "name of input elem file")
    ("cnnx", po::value<std::string>(), "name of input cnnx file")
    ("lon", po::value<std::string>(), "name of input lon file")
    ("output", po::value<std::string>(), "name of output file")
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

  if (!vm.count("pts")) 
  { 
    std::cout << "Input was not set.\n";
    return EXIT_FAILURE;
  } 

  if (!vm.count("output")) 
  { 
    std::cout << "Output was not set.\n";
    return EXIT_FAILURE;
  } 

  //Create vtk data structure and writer
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  create_points(vm["pts"].as<std::string>(),
                grid);

  if(vm.count("tris"))
  {
    create_tris(vm["tris"].as<std::string>(),
                grid);
  }

  if(vm.count("elem"))
  {
    create_elem(vm["elem"].as<std::string>(),
                grid);
  }

  if(vm.count("cnnx"))
  {
    create_lines(vm["cnnx"].as<std::string>(),
                 grid);
  }
 
  if(vm.count("lon"))
  {
    create_vectors(vm["lon"].as<std::string>(),
                   grid);
  }
  
  vtkXMLUnstructuredGridWriter* writer= vtkXMLUnstructuredGridWriter::New();
  writer->SetInput(grid);
  //writer->SetDataModeToAscii();
  writer->SetFileName(vm["output"].as<std::string>().c_str());
  writer->Write();
  writer->Delete();

}

//
//
//
void create_points(const std::string& file_name, 
                   vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  //parse file and simultaneously output to vtk
  std::ifstream pts_file(file_name.c_str());
  if (!pts_file)
  {
    std::cerr << "Failed to read pts file" << std::endl;
    abort();
  }

  //Read the pts file
  std::string line;
  std::getline(pts_file, line); //Skip over the first line

  //Create Vertices
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  unsigned vertex_index = 0;
  while(std::getline(pts_file, line))
  {
    double coords[3];
      
    char_delimiters_separator<char> sep(false, "", " ");
    tokenizer<> line_toks(line, sep);
      
    std::stringstream vertex_coords;

    tokenizer<>::iterator beg=line_toks.begin();
    std::istringstream iss(*beg);
    iss >> coords[0];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> coords[1];
    ++beg;

    std::istringstream iss3(*beg);
    iss3 >> coords[2];
      
    points->InsertPoint(vertex_index, coords);
    ++vertex_index;
  }

  grid->SetPoints(points);
}

//
//
//
void create_tris(const std::string& file_name,
                 vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  //Create tris
  std::ifstream tris_file(file_name.c_str());
  if (!tris_file)
  {
    std::cerr << "Failed to read tris file" << std::endl;
    abort();
  }

  //Read the pts file
  std::string line;
  std::getline(tris_file, line); //Skip over the first line

  //Create Vertices
  unsigned tri_index = 0;
  while(std::getline(tris_file, line))
  {
    std::vector<unsigned int> conn(3);
      
    char_delimiters_separator<char> sep(false, "", " ");
    tokenizer<> line_toks(line, sep);
      
    std::stringstream vertex_coords;

    tokenizer<>::iterator beg=line_toks.begin();
    std::istringstream iss(*beg);
    iss >> conn[0];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> conn[1];
    ++beg;

    std::istringstream iss3(*beg);
    iss3 >> conn[2];
    ++beg;
      
    vtkIdList *pts = vtkIdList::New();
    pts->SetNumberOfIds(3);

    for(unsigned int i=0;i<conn.size();++i)
    {
	  pts->SetId(i,conn[i]);
    } 
  
    grid->InsertNextCell(VTK_TRIANGLE,pts);
    pts->Delete();

    ++tri_index;
  }
}

//
//
//
void create_elem(const std::string& file_name,
                 vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  //Create elem
  std::ifstream elem_file(file_name.c_str());
  if (!elem_file)
  {
    std::cerr << "Failed to read elem file" << std::endl;
    abort();
  }

  //Read the pts file
  std::string line;
  std::getline(elem_file, line); //Skip over the first line

  //Create Vertices
  unsigned elem_index = 0;
  while(std::getline(elem_file, line))
  {
    std::vector<unsigned int> conn(4);
      
    char_delimiters_separator<char> sep(false, "", " ");
    tokenizer<> line_toks(line, sep);
      
    std::stringstream vertex_coords;

    tokenizer<>::iterator beg=line_toks.begin();
    ++beg; //Throw away Tt first entry
    std::istringstream iss(*beg);
    iss >> conn[0];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> conn[1];
    ++beg;

    std::istringstream iss3(*beg);
    iss3 >> conn[2];
    ++beg;
  
    std::istringstream iss4(*beg);
    iss4 >> conn[3];
    ++beg;
      
    vtkIdList *pts = vtkIdList::New();
    pts->SetNumberOfIds(4);

    for(unsigned int i=0;i<conn.size();++i)
    {
      //std::cout << i << " " << conn[i] << std::endl;
	  pts->SetId(i,conn[i]);
    //std::cout << pts->GetId(i) << std::endl;
    } 
  
    grid->InsertNextCell(VTK_TETRA,pts);
    pts->Delete();

    ++elem_index;
  }
}


/**
 * Reads a .cnnx file and adds the lines to the given grid;
 */
void create_lines(const std::string& file_name,
                  vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  //Create tris
  std::ifstream cnnx_file(file_name.c_str());
  if (!cnnx_file)
  {
    std::cerr << "Failed to read cnnx file" << std::endl;
    abort();
  }

  //Read the pts file
  std::string line;
  std::getline(cnnx_file, line); //Skip over the first line

  //Create Vertices
  unsigned cnnx_index = 0;
  while(std::getline(cnnx_file, line))
  {
    std::vector<unsigned int> conn(2);
      
    char_delimiters_separator<char> sep(false, "", " ");
    tokenizer<> line_toks(line, sep);
      
    std::stringstream vertex_coords;

    tokenizer<>::iterator beg=line_toks.begin();
    std::istringstream iss(*beg);
    iss >> conn[0];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> conn[1];
    ++beg;
      
    vtkIdList *pts = vtkIdList::New();
    pts->SetNumberOfIds(2);

    for(unsigned int i=0;i<conn.size();++i)
    {
	  pts->SetId(i,conn[i]);
    } 
  
    grid->InsertNextCell(VTK_LINE,pts);
    pts->Delete();

    ++cnnx_index;
  }
}


/**
 * Reads a .lon file and adds the vectors to the given elements;
 * Code assumes that the file has a first line to ignore.
 */
void create_vectors(const std::string& file_name,
                    vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  //Create tris
  std::ifstream lon_file(file_name.c_str());
  if (!lon_file)
  {
    std::cerr << "Failed to read lon file" << std::endl;
    abort();
  }
  
  // 
  
  //Read the pts file
  std::string line;
  // std::getline(lon_file, line); //Skip over the first line
  
  // quick and dirty way to check if there is a one-line header
  // if there is, ignore it
  // if (found!=string::npos)
  // 
  vtkCellData* cellData = grid->GetCellData();
  vtkSmartPointer< vtkFloatArray > vectors = vtkSmartPointer< vtkFloatArray >::New();
  vectors->SetName("lon");
  cellData->AddArray(vectors);
  
  //Create Vertices
  while(std::getline(lon_file, line))
  {
    double lon[3];
    
     // >> lon[0] >> lon[1] >> lon[2];
    
    char_delimiters_separator<char> sep(false, "", " ");
    tokenizer<> line_toks(line, sep);
      
    std::stringstream vertex_coords;

    tokenizer<>::iterator beg=line_toks.begin();
    std::istringstream iss(*beg);
    iss >> lon[0];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> lon[1];
    ++beg;
    
    std::istringstream iss3(*beg);
    iss3 >> lon[2];
    ++beg;
    
    vectors->SetNumberOfComponents(3);
    vectors->InsertNextTuple3(lon[0],lon[1],lon[2]);
    
    // vectors->InsertNextTuple(lon[0]);
    
    
    //       
    //     vtkIdList *pts = vtkIdList::New();
    //     pts->SetNumberOfIds(2);
    // 
    //     for(unsigned int i=0;i<lon.size();++i)
    //     {
    // pts->SetId(i,lon[i]);
    //     } 
    //   
    //     grid->InsertNextCell(VTK_LINE,pts);
    //     pts->Delete();
  }
}
