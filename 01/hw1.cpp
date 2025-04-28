#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/convex_hull_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

struct Face {
  unsigned int id;
  Kernel::Triangle_3 triangle;
  int material;
  std::set<unsigned int> neighbours;
  int block_id;
};

struct Block {
  std::vector<Kernel::Point_2> points;
  std::vector<Kernel::Point_2> convex_hull;
  Kernel::RT min, max;
};

typedef std::vector<Face>::iterator Faces_iterator;
typedef std::pair<Faces_iterator, Faces_iterator> Faces_intersection;
typedef CGAL::Box_intersection_d::Box_with_handle_d<Kernel::RT, 2, Faces_iterator> Box;

const std::string input_file = "/Users/ken/Downloads/10-282-562-obj/10-282-562-LoD12-3D.obj";
const std::string output_file = "/Users/ken/Downloads/delft.obj";

struct Box_intersector {
  std::back_insert_iterator<std::vector<Faces_intersection>> back_inserter;
  
  Box_intersector(const std::back_insert_iterator<std::vector<Faces_intersection>> &bi) : back_inserter(bi) { }
  
  void operator()(const Box &a, const Box &b) {
    *back_inserter++ = Faces_intersection(a.handle(), b.handle());
  }
};

int main(int argc, const char *argv[]) {
  
  std::vector<Kernel::Point_3> input_vertices;
  std::vector<Face> input_faces;
  std::vector<Box> boxes;
  std::vector<Block> blocks;
  
  double expansion_threshold = 1.0;
  
  // Read file
  std::ifstream input_stream;
  input_stream.open(input_file);
  if (input_stream.is_open()) {
    std::string line;
    int current_material = -1;
    
    // Parse line by line
    while (getline(input_stream, line)) {
//      std::cout << line << std::endl;
      
      std::istringstream line_stream(line);
      std::string line_type;
      line_stream >> line_type;
      
      // Vertex
      if (line_type == "v") {
        double x, y, z;
        line_stream >> x >> y >> z;
        input_vertices.emplace_back(Kernel::Point_3(x, y, z));
      }
      
      // Face
      if (line_type == "f") {
        input_faces.emplace_back();
        input_faces.back().id = input_faces.size()-1;
        input_faces.back().material = current_material;
        input_faces.back().block_id = -1;
        std::vector<Kernel::Point_3> face_vertices;
        int v;
        while (!line_stream.eof()) {
          line_stream >> v;
          face_vertices.emplace_back(input_vertices[v-1]);
        } if (face_vertices.size() == 3) {
          // to do...
        }
      }
      
      // Material
      if (line_type == "usemtl") {
        line_stream >> current_material;
      }
    }
  }
  
  // Print vertices
//  int current_vertex = 0;
//  for (auto const &vertex: input_vertices) {
//    std::cout << "Vertex " << current_vertex++ << ": " << "(" << vertex << ")" << std::endl;
//  }

  // Print faces
//  for (Faces_iterator current_face = input_faces.begin(); current_face != input_faces.end(); ++current_face) {
//    std::cout << "Face " << current_face->id << ":" << std::endl;
//    for (int i = 0; i < 3; ++i) {
//      std::cout << "\t(" << current_face->triangle[i] << ")" << std::endl;
//    }
//  }
  
  // Generate expanded bounding boxes
  for (Faces_iterator current_face = input_faces.begin(); current_face != input_faces.end(); ++current_face) {
    CGAL::Bbox_3 face_box = current_face->triangle.bbox();
    CGAL::Bbox_2 expanded_box;
    // to do...
    boxes.push_back(Box(expanded_box, current_face));
  }
  
  // Compute neighbours
  std::vector<Faces_intersection> intersections;
  Box_intersector box_intersector(std::back_inserter(intersections));
  CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), box_intersector);
  std::cout << intersections.size() << " intersections found" << std::endl;
  for (auto const &intersection: intersections) {
    // to do...
  }
  
  // Print face neighbours
//  for (auto const &face: input_faces) {
//    std::cout << "Face " << face.id << " neighbours:";
//    for (auto const &neighbour: face.neighbours) {
//      std::cout << " " << neighbour;
//    } std::cout << std::endl;
//  }
  
  // Assign blocks to faces
  int current_block = 0;
  std::list<unsigned int> faces_in_current_block;
  for (auto &face: input_faces) {
    if (face.block_id != -1) continue;
    face.block_id = current_block;
    faces_in_current_block.push_back(face.id);
    while (!faces_in_current_block.empty()) {
      for (auto const &neighbour: input_faces[faces_in_current_block.front()].neighbours) {
        if (input_faces[neighbour].block_id == -1) {
          input_faces[neighbour].block_id = current_block;
          faces_in_current_block.push_back(neighbour);
        }
      } faces_in_current_block.pop_front();
    } ++current_block;
  } std::cout << current_block << " blocks found" << std::endl;
  
  
  // Create blocks and put points in each
  // to do...
  
  // Compute convex hulls and min/max
  for (auto &block: blocks) {
    CGAL::convex_hull_2(block.points.begin(), block.points.end(), std::back_inserter(block.convex_hull));
    // to do...
  }
  
  // Create output vertices and faces
  std::vector<Kernel::Triangle_3> output_faces;
  // to do...
  
  // Write output
  std::ofstream output_stream;
  output_stream.open(output_file);
  // to do...
  output_stream.close();
  
  return 0;
}
