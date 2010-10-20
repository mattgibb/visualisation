/*
 * main.cpp
 */
#include <cmath>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGenericGeometryFilter.h"
#include "vtkCellArray.h"
#include "vtkConfigure.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"


enum vertex_x_t {vertex_x = 1000};
enum vertex_y_t {vertex_y = 1001};
enum vertex_z_t {vertex_z = 1002};

namespace boost
{
	BOOST_INSTALL_PROPERTY(vertex, x);
	BOOST_INSTALL_PROPERTY(vertex, y);
	BOOST_INSTALL_PROPERTY(vertex, z);
}

using namespace boost;

// template <typename DistanceMap, typename EdgeMap>
// class my_distance_recorder : public default_bfs_visitor
// {
// public:
//   my_distance_recorder(DistanceMap dist, EdgeMap edge_map) : d(dist), _edge_map(edge_map) { }

//   template <typename Edge, typename Graph>
//   void tree_edge(Edge e, const Graph& g) const
//   {
//     typename graph_traits<Graph>::vertex_descriptor
//       u = source(e, g), v = target(e, g);
//     d[v] = d[u] + _edge_map[e];
//   }
// private:
//   DistanceMap d;
//   EdgeMap _edge_map;
// };

// Convenience function
// template <typename DistanceMap, typename EdgeMap>
// my_distance_recorder<DistanceMap, EdgeMap> record_distance (DistanceMap d, EdgeMap e)
// {
//   return my_distance_recorder<DistanceMap, EdgeMap>(d, e); 
// }

// auxiliary types
struct location
{
	float x, y, z; // coordinates
};

//Boost graph typedefs
typedef adjacency_list<vecS, vecS, undirectedS,
 					  property<vertex_x_t, float,
 					  property<vertex_y_t, float,
 					  property<vertex_z_t, float > > >,
 					  property<edge_weight_t, float > > Graph;

typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef graph_traits < Graph >::edge_descriptor Edge;
typedef std::map < unsigned, Vertex > IndexVertexMap;


namespace po = boost::program_options;

//Prototypes

/**
 * Write the resulting graph out to vtk
 */
void write_to_vtk(Graph& g, 
                  IndexVertexMap& vertices,
                  property_map<Graph, vertex_index_t>::type& vertex_index_map,
                  property_map<Graph, vertex_x_t>::type& vertex_x_map,
                  property_map<Graph, vertex_y_t>::type& vertex_y_map,
                  property_map<Graph, vertex_z_t>::type& vertex_z_map,
                  std::vector<int>& distance,
                  const std::string& pkj_input_root,
                  const std::string& surface_input_root,
                  const unsigned& vertex_offset,
                  const unsigned& edge_offset);

/**
 * Parse a set of files specified by pkj_input_root to create the pkj portion of the graph
 */
void parse_pkj(Graph& g, 
               IndexVertexMap& vertices,
               property_map<Graph, vertex_index_t>::type& vertex_index_map,
               property_map<Graph, vertex_x_t>::type& vertex_x_map,
               property_map<Graph, vertex_y_t>::type& vertex_y_map,
               property_map<Graph, vertex_z_t>::type& vertex_z_map,
               property_map<Graph, edge_weight_t>::type& edge_length_map,
               const std::string& input_root);

/**
 * Parse the pmj file
 */
void parse_pmj(Graph& g, 
               IndexVertexMap& vertices,
               property_map<Graph, vertex_index_t>::type& vertex_index_map,
               property_map<Graph, vertex_x_t>::type& vertex_x_map,
               property_map<Graph, vertex_y_t>::type& vertex_y_map,
               property_map<Graph, vertex_z_t>::type& vertex_z_map,
               property_map<Graph, edge_weight_t>::type& edge_length_map,
               unsigned offset,
               const std::string& input_root);

/**
 * Parse an unstructured grid file representing the surface
 */
void parse_surface(Graph& g, 
                   IndexVertexMap& vertices,
                   property_map<Graph, vertex_index_t>::type& vertex_index_map,
                   property_map<Graph, vertex_x_t>::type& vertex_x_map,
                   property_map<Graph, vertex_y_t>::type& vertex_y_map,
                   property_map<Graph, vertex_z_t>::type& vertex_z_map,
                   property_map<Graph, edge_weight_t>::type& edge_length_map,
                   const std::string& input_root);

int main(int argc, char *argv[])
{
  //Process command line 
  //set up possible options
  po::options_description desc("Command line arguments");
  desc.add_options()
    ("help", "produce help message")
    ("pkj", po::value<std::string>(), "Root of free-running pkj input file")
    ("surface", po::value<std::string>(), "Endocardial surface unstructured grid file")
    ("pmj", po::value<std::string>(), "Pmjs input file")
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

  if (!vm.count("pkj")) 
  { 
    std::cout << "pkj input was not set.\n";
    return EXIT_FAILURE;
  } 

  if (!vm.count("surface")) 
  { 
    std::cout << "surface input was not set.\n";
    return EXIT_FAILURE;
  } 

  if (!vm.count("pmj")) 
  { 
    std::cout << "pmj input was not set.\n";
    return EXIT_FAILURE;
  } 

  std::string pkj_input_root = vm["pkj"].as<std::string>();
  std::cout << "Pkj input file was set to: "
            << pkj_input_root << ".\n";

  std::string surface_input_root = vm["surface"].as<std::string>();
  std::cout << "Surface input file was set to: "
            << surface_input_root << ".\n";

  std::string pmj_input_root = vm["pmj"].as<std::string>();
  std::cout << "Pmj input file was set to: "
            << pmj_input_root << ".\n";


  //Create the graph & associated types
  Graph g;
  IndexVertexMap vertices;

  property_map<Graph, vertex_index_t>::type vertex_index_map = get(vertex_index, g);
  property_map<Graph, vertex_x_t>::type vertex_x_map = get(vertex_x, g);
  property_map<Graph, vertex_y_t>::type vertex_y_map = get(vertex_y, g);
  property_map<Graph, vertex_z_t>::type vertex_z_map = get(vertex_z, g);

  property_map<Graph, edge_weight_t>::type edge_length_map = get(edge_weight, g);

  parse_pkj(g, 
            vertices, 
            vertex_index_map, 
            vertex_x_map, 
            vertex_y_map, 
            vertex_z_map, 
            edge_length_map,
            pkj_input_root);

  unsigned vertex_offset = num_vertices(g);
  unsigned edge_offset = num_edges(g);

  parse_surface(g, 
                vertices, 
                vertex_index_map, 
                vertex_x_map, 
                vertex_y_map, 
                vertex_z_map, 
                edge_length_map,
                surface_input_root);


  parse_pmj(g, 
            vertices, 
            vertex_index_map, 
            vertex_x_map, 
            vertex_y_map, 
            vertex_z_map, 
            edge_length_map,
            vertex_offset,
            pmj_input_root);

  //Record the distances in a vector
  std::vector<int> distance(num_vertices(g));
  std::vector<Vertex> predecessor(num_vertices(g));

  Vertex src = vertices[55908 + vertex_offset];
  distance[src] = 0;

  dijkstra_shortest_paths(g, src, predecessor_map(&predecessor[0]).distance_map(&distance[0]));
  // breadth_first_search(g, src, visitor(record_distance(&distance[0], edge_length_map)));

  std::cout << "search done" << std::endl;

  //Now loop through the graph and output the results...
  write_to_vtk(g,
               vertices, 
               vertex_index_map, 
               vertex_x_map, 
               vertex_y_map, 
               vertex_z_map, 
               distance,
               pkj_input_root,
               surface_input_root,
               vertex_offset,
               edge_offset);

}

//
//
//
void write_to_vtk(Graph& g, 
                  IndexVertexMap& index_vertex_map,
                  property_map<Graph, vertex_index_t>::type& vertex_index_map,
                  property_map<Graph, vertex_x_t>::type& vertex_x_map,
                  property_map<Graph, vertex_y_t>::type& vertex_y_map,
                  property_map<Graph, vertex_z_t>::type& vertex_z_map,
                  std::vector<int>& distance,
                  const std::string& pkj_input_root,
                  const std::string& surface_input_root,
                  const unsigned& offset,
                  const unsigned& edge_offset)
{
  vtkUnstructuredGrid* _vtk_grid = vtkUnstructuredGrid::New();
  vtkXMLPUnstructuredGridWriter* writer= vtkXMLPUnstructuredGridWriter::New();
 
  //Output the nodes of the Purkinje network
  vtkPoints* points = vtkPoints::New();
  vtkFloatArray *data = vtkFloatArray::New();
  data->SetName("distance");

  //Loop over vertices...
  graph_traits<Graph>::vertex_iterator i, end;

  unsigned vertex_index = 0;

  for (boost::tie(i, end) = vertices(g); i != end && vertex_index < offset; ++i)
  {
    //Set nodes
    double coords[3];
    coords[0] = vertex_x_map[*i];
    coords[1] = vertex_y_map[*i];
    coords[2] = vertex_z_map[*i];

    points->InsertPoint(vertex_index_map[*i], coords);

    //Output distance data
    data->InsertValue(vertex_index_map[*i], distance[vertex_index_map[*i]]);

    ++vertex_index;
  }
  
  _vtk_grid->SetPoints(points);
  _vtk_grid->GetPointData()->AddArray(data);

  //Output connections
  graph_traits<Graph>::edge_iterator j, endj;

  unsigned edge_index = 0;
  for (boost::tie(j, endj) = edges(g); j != endj && edge_index < edge_offset; ++j)
  {
    vtkIdList *pts = vtkIdList::New();
    pts->SetNumberOfIds(2);
    pts->SetId(0, source(*j, g));
    pts->SetId(1, target(*j, g));

    _vtk_grid->InsertNextCell(VTK_LINE,pts);

    ++edge_index;

    pts->Delete();
  }


  writer->SetInput(_vtk_grid);
  writer->SetFileName((pkj_input_root + ".vtk").c_str());
  writer->Write();

  data->Delete();
  points->Delete();
  _vtk_grid->Delete();

  //Now load up the original surface file 
  vtkXMLUnstructuredGridReader* surface_reader = vtkXMLUnstructuredGridReader::New();
  surface_reader->SetFileName((surface_input_root + ".vtu").c_str());
  surface_reader->Update();
  vtkUnstructuredGrid* surface_grid = surface_reader->GetOutput();

  vtkFloatArray *surface_data = vtkFloatArray::New();
  surface_data->SetName("distance");

  for (; i != end; ++i)
  {
    //Output distance data
    //std::cout << vertex_index_map[*i] - offset << std::endl;
    surface_data->InsertValue(vertex_index_map[*i] - offset, distance[vertex_index_map[*i]]);
  }

  
  //  std::cout << offset << std::endl;
  //std::cout << surface_data->GetSize() << std::endl;
  //std::cout << surface_grid->GetPoints()->GetData()->GetSize()/3 << std::endl;

  surface_grid->GetPointData()->AddArray(surface_data);

  writer->SetInput(surface_grid);
  writer->SetFileName((surface_input_root + "_out.vtk").c_str());
  writer->Write();

  surface_data->Delete();
  surface_reader->Delete();
  writer->Delete();
}

//
//
//
void parse_pkj(Graph& g, 
               IndexVertexMap& vertices,
               property_map<Graph, vertex_index_t>::type& vertex_index_map,
               property_map<Graph, vertex_x_t>::type& vertex_x_map,
               property_map<Graph, vertex_y_t>::type& vertex_y_map,
               property_map<Graph, vertex_z_t>::type& vertex_z_map,
               property_map<Graph, edge_weight_t>::type& edge_length_map,
               const std::string& input_root)
{
  std::ifstream pts_file((input_root + ".pts").c_str());
  if (!pts_file)
  {
	std::cerr << "Failed to read pts file" << std::endl;
	abort();
  }

  //Read the pts file
  std::string line;
  std::getline(pts_file, line); //Skip over the first line

  //Create Vertices
  while(std::getline(pts_file, line))
  {
	char_delimiters_separator<char> sep(false, "", " ");
	tokenizer<> line_toks(line, sep);

    std::stringstream vertex_coords;

    Vertex u;
    u = add_vertex(g);

    tokenizer<>::iterator beg=line_toks.begin();
    std::istringstream iss(*beg);
    iss >> vertex_x_map[u];
    ++beg;

    std::istringstream iss2(*beg);
    iss2 >> vertex_y_map[u];
    ++beg;

    std::istringstream iss3(*beg);
    iss3 >> vertex_z_map[u];

    IndexVertexMap::iterator pos;
    bool inserted;

    tie(pos, inserted) = vertices.insert(std::make_pair(vertex_index_map[u], u));

    //std::cout << vertex_x_map[u] << " " << vertex_y_map[u] << " " << vertex_z_map[u] << "\n";
  }

  pts_file.close();

  //Read the cnnx file

	// //Create Edges
  std::ifstream cnnx_file((input_root + ".cnnx").c_str());
  if (!cnnx_file)
  {
	std::cerr << "Failed to read cnnx file" << std::endl;
	abort();
  }

  std::getline(cnnx_file, line); //Skip over the first line

  unsigned i = 0;

  //Create Edges
  while(std::getline(cnnx_file, line))
  {
     char_delimiters_separator<char> sep(false, "", " ");
     tokenizer<> line_toks(line, sep);

     unsigned pt0;
     unsigned pt1;

     tokenizer<>::iterator beg=line_toks.begin();
     std::istringstream iss(*beg);
     iss >> pt0;
     ++beg;

     std::istringstream iss2(*beg);
     iss2 >> pt1;
     ++beg;

	 ++i;

     Vertex p0 = vertices[pt0];
     Vertex p1 = vertices[pt1];
 
     bool inserted;
 
     graph_traits<Graph>::edge_descriptor e0;
     tie(e0, inserted) = add_edge(p0, p1, g);
     edge_length_map[e0] = std::sqrt(std::pow(vertex_x_map[p0] - vertex_x_map[p1], 2) +
							   		   std::pow(vertex_y_map[p0] - vertex_y_map[p1], 2) +
							   		   std::pow(vertex_z_map[p0] - vertex_z_map[p1], 2));

     //if p0 & p1 are start nodes then set edge map length to 0..., else use Euclidean distance
     //could just check location and use a region based approach...
  }
  
  cnnx_file.close();


  //Read the second pts file...
  //Load all points as vertices and record offset...

  //Load the face file as connections, mung connections that are within the start zone/range...

  //How will we add pmjs? This should be done offline...
  
  return;
}


void parse_surface(Graph& g, 
                   IndexVertexMap& vertices,
                   property_map<Graph, vertex_index_t>::type& vertex_index_map,
                   property_map<Graph, vertex_x_t>::type& vertex_x_map,
                   property_map<Graph, vertex_y_t>::type& vertex_y_map,
                   property_map<Graph, vertex_z_t>::type& vertex_z_map,
                   property_map<Graph, edge_weight_t>::type& edge_length_map,
                   const std::string& input_root)
{
  //Septum bounding Coordinates
  double xmin = 1.01e4;
  double xmax = 1.4e4;
  double ymin = 8.82e3;
  double ymax = 1.59e4;
  double zmin = 7.13e3;
  double zmax = 2.3e4;


  vtkXMLUnstructuredGridReader* surface_reader = vtkXMLUnstructuredGridReader::New();
  surface_reader->SetFileName((input_root + ".vtu").c_str());
  surface_reader->Update();
  vtkUnstructuredGrid* surface_grid = surface_reader->GetOutput();
  vtkPoints* surface_points = surface_grid->GetPoints(); 
  vtkCellArray* surface_cells = surface_grid->GetCells();

  //Work out current offset...
  unsigned vertex_offset = num_vertices(g);

  //Loop over points file and create vertices
  for (vtkIdType i = 0; i < surface_points->GetNumberOfPoints(); ++i)
  {
    double coords[3];
    surface_points->GetPoint(i, coords);

    Vertex u;
    u = add_vertex(g);

    vertex_x_map[u] = coords[0];
    vertex_y_map[u] = coords[1];
    vertex_z_map[u] = coords[2];

    IndexVertexMap::iterator pos;
    bool inserted;

    tie(pos, inserted) = vertices.insert(std::make_pair(vertex_index_map[u], u));
  }

  //Loop over the cells and create the connections
  vtkIdType npts;
  vtkIdType* pts;
  for(surface_cells->InitTraversal(); surface_cells->GetNextCell(npts, pts);)
  {
     Vertex p0 = vertices[pts[0] + vertex_offset];
     Vertex p1 = vertices[pts[1] + vertex_offset];
     Vertex p2 = vertices[pts[2] + vertex_offset];
 
     bool p0_in_septum = vertex_x_map[p0] > xmin && vertex_x_map[p0] < xmax &&
                         vertex_y_map[p0] > ymin && vertex_y_map[p0] < ymax &&
                         vertex_z_map[p0] > zmin && vertex_z_map[p0] < zmax;
     bool p1_in_septum = vertex_x_map[p1] > xmin && vertex_x_map[p1] < xmax &&
                         vertex_y_map[p1] > ymin && vertex_y_map[p1] < ymax &&
                         vertex_z_map[p1] > zmin && vertex_z_map[p1] < zmax;
     bool p2_in_septum = vertex_x_map[p2] > xmin && vertex_x_map[p2] < xmax &&
                         vertex_y_map[p2] > ymin && vertex_y_map[p2] < ymax &&
                         vertex_z_map[p2] > zmin && vertex_z_map[p2] < zmax;

     bool inserted;
 
     graph_traits<Graph>::edge_descriptor e0;
     tie(e0, inserted) = add_edge(p0, p1, g);
     
     if(p0_in_septum && p1_in_septum)
     {
       edge_length_map[e0] = 0.0;
     }
     else
     {
       edge_length_map[e0] = std::sqrt(std::pow(vertex_x_map[p0] - vertex_x_map[p1], 2) +
							   		   std::pow(vertex_y_map[p0] - vertex_y_map[p1], 2) +
							   		   std::pow(vertex_z_map[p0] - vertex_z_map[p1], 2));
     }

     graph_traits<Graph>::edge_descriptor e1;
     tie(e1, inserted) = add_edge(p0, p2, g);

     if(p0_in_septum && p2_in_septum)
     {
       edge_length_map[e1] = 0.0;
     }
     else
     {
       edge_length_map[e1] = std::sqrt(std::pow(vertex_x_map[p0] - vertex_x_map[p2], 2) +
							   		   std::pow(vertex_y_map[p0] - vertex_y_map[p2], 2) +
							   		   std::pow(vertex_z_map[p0] - vertex_z_map[p2], 2));
     }

     graph_traits<Graph>::edge_descriptor e2;
     tie(e2, inserted) = add_edge(p1, p2, g);
     
     if(p1_in_septum && p2_in_septum)
     {
       edge_length_map[e2] = 0.0;
     }
     else
     {
       edge_length_map[e2] = std::sqrt(std::pow(vertex_x_map[p1] - vertex_x_map[p2], 2) +
							   		   std::pow(vertex_y_map[p1] - vertex_y_map[p2], 2) +
							   		   std::pow(vertex_z_map[p1] - vertex_z_map[p2], 2));
     }

  }
  
  surface_reader->Delete();
}

void parse_pmj(Graph& g, 
               IndexVertexMap& vertices,
               property_map<Graph, vertex_index_t>::type& vertex_index_map,
               property_map<Graph, vertex_x_t>::type& vertex_x_map,
               property_map<Graph, vertex_y_t>::type& vertex_y_map,
               property_map<Graph, vertex_z_t>::type& vertex_z_map,
               property_map<Graph, edge_weight_t>::type& edge_length_map,
               unsigned offset,
               const std::string& input_root)
{

	// //Create Edges
  std::ifstream pmj_file(input_root.c_str());
  if (!pmj_file)
  {
	std::cerr << "Failed to read pmj file" << std::endl;
	abort();
  }

  //  std::getline(cnnx_file, line); //Skip over the first line

  unsigned i = 0;

  //Create Edges
  std::string line;
  while(std::getline(pmj_file, line))
  {
     char_delimiters_separator<char> sep(false, "", " ");
     tokenizer<> line_toks(line, sep);

     unsigned pt0;
     unsigned pt1;

     tokenizer<>::iterator beg=line_toks.begin();
     std::istringstream iss(*beg);
     iss >> pt0;
     ++beg;

     std::istringstream iss2(*beg);
     iss2 >> pt1;
     ++beg;

	 ++i;

     pt1 = pt1 + offset;
     //std::cout << "pt0: " << pt0 << " pt1: " << pt1 << std::endl;

     Vertex p0 = vertices[pt0];
     Vertex p1 = vertices[pt1];
     //std::cout << vertices.size() << std::endl;

     //std::cout << "pkj: " << vertex_x_map[p0] << " "
     //          << vertex_y_map[p0] << " " 
     //          << vertex_z_map[p0] << " " 
     //          << "sur: "<< vertex_x_map[p1] << " "
     //          << vertex_y_map[p1] << " " 
     //          << vertex_z_map[p1] << std::endl;
     
     bool inserted;
 
     graph_traits<Graph>::edge_descriptor e0;
     tie(e0, inserted) = add_edge(p0, p1, g);
     edge_length_map[e0] = std::sqrt(std::pow(vertex_x_map[p0] - vertex_x_map[p1], 2) +
     						   		   std::pow(vertex_y_map[p0] - vertex_y_map[p1], 2) +
     						   		   std::pow(vertex_z_map[p0] - vertex_z_map[p1], 2));
  }

  std::cout << offset << std::endl;
}
