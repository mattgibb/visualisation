//
// Generates bundle branches from a set of control point and a set of "pmjs"
//
#include <sstream>
#include <cmath>

#include <boost/program_options.hpp>

#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkXMLPPolyDataWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkConfigure.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkSmartPointer.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkMath.h"

using namespace boost;

namespace po = boost::program_options;


/**
 * Reads a pmjs file and return (filtered) end points
 */

void read_pmjs(po::variables_map& vm,
               const std::string& prefix,
               std::vector<vtkIdType>& pmjs,
               vtkSmartPointer<vtkUnstructuredGrid>& surface_grid);

/**
 * Writes a shortest geodesic path to an unstructured grid.
 */
void create_shortest_path(vtkSmartPointer<vtkGeometryFilter> geometry,
                          vtkSmartPointer<vtkDijkstraGraphGeodesicPath> geodesic,
                          vtkSmartPointer<vtkFloatArray> distance,
                          vtkSmartPointer<vtkUnstructuredGrid> pkj_grid,
                          const vtkIdType& branch_id);

/**
 * Writes bundle branches to an unstructured grid
 */
void create_bundle_branches(vtkSmartPointer<vtkUnstructuredGrid> surface_grid,
                            vtkSmartPointer<vtkUnstructuredGrid> pkj_grid,
                            vtkSmartPointer<vtkFloatArray> distance,
                            const vtkIdType& his_surface_id,
                            const vtkIdType& lv_branch_surface_id,
                            const vtkIdType& rv_branch_surface_id,
                            vtkIdType& lv_branch_id,
                            vtkIdType& rv_branch_id);

int main(int argc, char *argv[])
{
  //Process command line 
  //set up possible options
  po::options_description desc("Command line arguments");
  desc.add_options()
    ("help", "produce help message")
    ("surface", po::value<std::string>(), "Endocardial surface unstructured grid file")
    ("output", po::value<std::string>(), "Output file name")
    ("lv_pmj", po::value<std::string>(), "Pmjs input file for the left ventricle")
    ("rv_pmj", po::value<std::string>(), "Pmjs input file for the right ventricle")
    ("his", po::value<unsigned>(), "node index position for the bundle of his")
    ("lv_branch", po::value<unsigned>(), "node index position for the lv bundle branch")
    ("rv_branch", po::value<unsigned>(), "node index position for the rv bundle branch")
    ("lv_septum_xmin", po::value<double>(), "Definition of the lv septal region")
    ("lv_septum_xmax", po::value<double>(), "Definition of the lv septal region")
    ("lv_septum_ymin", po::value<double>(), "Definition of the lv septal region")
    ("lv_septum_ymax", po::value<double>(), "Definition of the lv septal region")
    ("lv_septum_zmin", po::value<double>(), "Definition of the lv septal region")
    ("lv_septum_zmax", po::value<double>(), "Definition of the lv septal region")
    ("rv_septum_xmin", po::value<double>(), "Definition of the rv septal region")
    ("rv_septum_xmax", po::value<double>(), "Definition of the rv septal region")
    ("rv_septum_ymin", po::value<double>(), "Definition of the rv septal region")
    ("rv_septum_ymax", po::value<double>(), "Definition of the rv septal region")
    ("rv_septum_zmin", po::value<double>(), "Definition of the rv septal region")
    ("rv_septum_zmax", po::value<double>(), "Definition of the rv septal region")
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

  std::string surface_file = vm["surface"].as<std::string>();

  //read surface file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> surface_reader = 
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  surface_reader->SetFileName(surface_file.c_str());
  surface_reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> surface_grid = surface_reader->GetOutput();
  vtkSmartPointer<vtkPoints> surface_points = surface_grid->GetPoints();

  //
  //Create the bundle of his & the bundle branches in a separate unstructured grid file
  //
  vtkSmartPointer<vtkUnstructuredGrid> pkj_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  pkj_grid->SetPoints(vtkSmartPointer<vtkPoints>::New());
  vtkSmartPointer<vtkFloatArray> distance = vtkSmartPointer<vtkFloatArray>::New();
  distance->SetName("distance");
  pkj_grid->GetPointData()->AddArray(distance);
  
  vtkIdType lv_branch_id;
  vtkIdType rv_branch_id;
  create_bundle_branches(surface_grid, 
                         pkj_grid,
                         distance,
                         vm["his"].as<unsigned>(),
                         vm["lv_branch"].as<unsigned>(),
                         vm["rv_branch"].as<unsigned>(),
                         lv_branch_id,
                         rv_branch_id);
  
  //Convert the unstructured grid geometry to a poly data geometry for use with geodesic
  vtkSmartPointer<vtkGeometryFilter> geometry = vtkSmartPointer<vtkGeometryFilter>::New();
  geometry->SetInput(surface_grid);
  geometry->Update();

  //Set up a geodesic path finder
  vtkSmartPointer<vtkDijkstraGraphGeodesicPath> geodesic =
     vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
  
  geodesic->SetInput(geometry->GetOutput());

  //Do lv branches, get a list of pmjs to connect up
  {
    std::vector<vtkIdType> pmjs;
    read_pmjs(vm, "lv", pmjs, surface_grid);

    for (std::vector<vtkIdType>::iterator iter = pmjs.begin();
         iter != pmjs.end();
         ++iter)
    {
      geodesic->SetStartVertex(vm["lv_branch"].as<unsigned>());
      geodesic->SetEndVertex(*iter);
      geodesic->SetStopWhenEndReached(true);
      geodesic->Update();

      create_shortest_path(geometry,
                           geodesic,
                           distance,
                           pkj_grid,
                           lv_branch_id);
    }
  }

  //Do rv branches
  {
    std::vector<vtkIdType> pmjs;
    read_pmjs(vm, "rv", pmjs, surface_grid);

    for (std::vector<vtkIdType>::iterator iter = pmjs.begin();
         iter != pmjs.end();
         ++iter)
    {
      geodesic->SetStartVertex(vm["rv_branch"].as<unsigned>());
      geodesic->SetEndVertex(*iter);
      geodesic->SetStopWhenEndReached(true);
      geodesic->Update();

      create_shortest_path(geometry,
                           geodesic,
                           distance,
                           pkj_grid,
                           rv_branch_id);
    }
  }

  //Output the generated fibres
  vtkXMLPUnstructuredGridWriter* writer= vtkXMLPUnstructuredGridWriter::New();
  writer->SetInput(pkj_grid);
  writer->SetFileName(vm["output"].as<std::string>().c_str());
  writer->Write();
  writer->Delete();
}

void read_pmjs(po::variables_map& vm,
               const std::string& prefix,
               std::vector<vtkIdType>& pmjs,
               vtkSmartPointer<vtkUnstructuredGrid>& surface_grid)
{
  //load the relevant file
  std::ifstream pmjs_file(vm[prefix + "_pmj"].as<std::string>().c_str());
  if (!pmjs_file)
  {
	std::cerr << "Failed to read .pmj file" << std::endl;
	abort();
  }

  //Read the file
  std::string line;

  while(std::getline(pmjs_file, line))
  {
	char_delimiters_separator<char> sep(false, "", " ");
	tokenizer<> line_toks(line, sep);
    
    vtkIdType pkj_pt_id;

    tokenizer<>::iterator beg=line_toks.begin();
    ++beg; //skip the first entry
    std::istringstream iss(*beg);
    iss >> pkj_pt_id; 

    
    //check validity
    double location[3];
    surface_grid->GetPoints()->GetPoint(pkj_pt_id, location);

    std::cout << pkj_pt_id << " " << location[0] << " " << location[1] << " " << location[2] << std::endl;

    // get coords, compare with given coords 
    if (location[0] > vm[prefix + "_septum_xmin"].as<double>() && 
        location[0] < vm[prefix + "_septum_xmax"].as<double>() &&
        location[1] > vm[prefix + "_septum_ymin"].as<double>() && 
        location[1] < vm[prefix + "_septum_ymax"].as<double>() &&
        location[2] > vm[prefix + "_septum_zmin"].as<double>() && 
        location[2] < vm[prefix + "_septum_zmax"].as<double>())
    {
      // if so, add pkj_pt_id to pmjs
      pmjs.push_back(pkj_pt_id);
      std::cout << "in!" << std::endl;
    }

  }
}


void create_shortest_path(vtkSmartPointer<vtkGeometryFilter> geometry,
                          vtkSmartPointer<vtkDijkstraGraphGeodesicPath> geodesic,
                          vtkSmartPointer<vtkFloatArray> distance,
                          vtkSmartPointer<vtkUnstructuredGrid> pkj_grid,
                          const vtkIdType& branch_id)
{
  vtkIdList* path_ids = geodesic->GetIdList();
  
  //Check if path failed (-1 in second entry)
  if(path_ids->GetId(1) == -1)
  {
    std::cout << "Path failed!" << std::endl;
    return;
  }

  //vtk returns the path from the target to the source, so we must build it in reverse
  double dist = distance->GetValue(branch_id);
  
  //Build the first branch
  double temp_point[3];
  geometry->GetOutput()->GetPoints()->GetPoint(path_ids->GetId(path_ids->GetNumberOfIds()-2), temp_point);  
  dist += std::sqrt(vtkMath::Distance2BetweenPoints(pkj_grid->GetPoint(branch_id), temp_point));

  vtkIdType start = pkj_grid->GetPoints()->InsertNextPoint(temp_point);
  distance->InsertNextValue(dist);
  
  //need to handle final line differently, as it connects to the bundle branch start
  vtkIdType ptIds[2];
  ptIds[0] = branch_id;
  ptIds[1] = start;
  pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);

  //
  // This doesn't quite output the whole line, for some reason i >= 0 lets i become < 0!
  //
  for (unsigned i = path_ids->GetNumberOfIds() - 3; i > 0; --i)
  { 
    //std::cout << i << " " << std::endl;
    //std::cout << path_ids->GetId(i+1) << std::endl;

    //create point + distance
    std::cout << "i: " << i << std::endl;
    geometry->GetOutput()->GetPoints()->GetPoint(path_ids->GetId(i), temp_point);  
    vtkIdType end = pkj_grid->GetPoints()->InsertNextPoint(temp_point);

    std::cout << i << " " << path_ids->GetId(i+1) << std::endl;
    std::cout << temp_point[0] << " " << temp_point[1] << " " << temp_point[2] << std::endl;
    double previous_point[3];
    geometry->GetOutput()->GetPoints()->GetPoint(path_ids->GetId(i+1), previous_point);
    dist += std::sqrt(
      vtkMath::Distance2BetweenPoints(previous_point, temp_point));

    std::cout << "dist:" << dist << std::endl;
    distance->InsertNextValue(dist);
    
    ptIds[0] = start;
    ptIds[1] = end;
    pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);

    start = end;
  }

  /*  double temp_point[3];
  geometry->GetOutput()->GetPoints()->GetPoint(path_ids->GetId(0), temp_point);

  vtkIdType start = pkj_grid->GetPoints()->InsertNextPoint(temp_point);

  for (unsigned i = 1; i < path_ids->GetNumberOfIds()-1; ++i)
  {
    geometry->GetOutput()->GetPoints()->GetPoint(path_ids->GetId(i), temp_point);
    vtkIdType end = pkj_grid->GetPoints()->InsertNextPoint(temp_point);

    vtkIdType ptIds[2];
    ptIds[0] = start;
    ptIds[1] = end;
    pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);

    start = end;
  }

  //need to handle final line differently, as it connects to the bundle branch start
  vtkIdType ptIds[2];
  ptIds[0] = start;
  ptIds[1] = branch_id;
  pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);*/
}


void create_bundle_branches(vtkSmartPointer<vtkUnstructuredGrid> surface_grid,
                            vtkSmartPointer<vtkUnstructuredGrid> pkj_grid,
                            vtkSmartPointer<vtkFloatArray> distance,
                            const vtkIdType& his_surface_id,
                            const vtkIdType& lv_branch_surface_id,
                            const vtkIdType& rv_branch_surface_id,
                            vtkIdType& lv_branch_id,
                            vtkIdType& rv_branch_id)
{
  //Get the coords of the points to create
  double his_coords[3], lv_branch_coords[3], rv_branch_coords[3], branch_coords[3];
  surface_grid->GetPoints()->GetPoint(his_surface_id, his_coords);
  surface_grid->GetPoints()->GetPoint(lv_branch_surface_id, lv_branch_coords);  
  surface_grid->GetPoints()->GetPoint(rv_branch_surface_id, rv_branch_coords);

  //Calculate the branch point location
  branch_coords[0] = (lv_branch_coords[0] + rv_branch_coords[0])/2;
  branch_coords[1] = (lv_branch_coords[1] + rv_branch_coords[1])/2;
  branch_coords[2] = (lv_branch_coords[2] + rv_branch_coords[2])/2;

  //create the points
  vtkIdType his_id = pkj_grid->GetPoints()->InsertNextPoint(his_coords);
  distance->InsertNextValue(0);

  double branch_dist = std::sqrt(vtkMath::Distance2BetweenPoints(his_coords,branch_coords));
  vtkIdType branch_id = pkj_grid->GetPoints()->InsertNextPoint(branch_coords);
  distance->InsertNextValue(branch_dist);

  lv_branch_id = pkj_grid->GetPoints()->InsertNextPoint(lv_branch_coords);
  distance->InsertNextValue(branch_dist +
                            std::sqrt(vtkMath::Distance2BetweenPoints(branch_coords, lv_branch_coords)));

  rv_branch_id = pkj_grid->GetPoints()->InsertNextPoint(rv_branch_coords);
  distance->InsertNextValue(branch_dist + 
                            std::sqrt(vtkMath::Distance2BetweenPoints(branch_coords,rv_branch_coords)));

  //create line segments
  vtkIdType ptIds[2];
  ptIds[0] = his_id;
  ptIds[1] = branch_id;
  pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);

  ptIds[0] = branch_id;
  ptIds[1] = lv_branch_id;
  pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);

  ptIds[0] = branch_id;
  ptIds[1] = rv_branch_id;
  pkj_grid->InsertNextCell(VTK_LINE, 2, ptIds);
}
