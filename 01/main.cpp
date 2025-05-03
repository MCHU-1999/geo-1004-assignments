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
  std::set<unsigned int> neighbours;  // std::set only stores unique elements
  int block_id;
};

struct Block {
  std::vector<Kernel::Point_2> points;
  std::vector<Kernel::Point_2> convex_hull;
  std::vector<Kernel::Point_3> convex_hull_3;
  Kernel::RT min=std::numeric_limits<double>::max(), max=std::numeric_limits<double>::lowest();
  // RT: Ring Type. Here they stands for z-min/max. 
  // RT is basically just number, but don't have to worry about what exact number type.
};

typedef std::vector<Face>::iterator Faces_iterator;
typedef std::pair<Faces_iterator, Faces_iterator> Faces_intersection;
typedef CGAL::Box_intersection_d::Box_with_handle_d<Kernel::RT, 2, Faces_iterator> Box;

// A struct with a callback function after an pair of intersecting boxes are identified.
struct Box_intersector {
  std::back_insert_iterator<std::vector<Faces_intersection>> back_inserter;
  
  Box_intersector(const std::back_insert_iterator<std::vector<Faces_intersection>> &bi) : back_inserter(bi) { }
  
  void operator()(const Box &a, const Box &b) {
    *back_inserter++ = Faces_intersection(a.handle(), b.handle());
  }
};

std::pair<double, double> triangle_z_minmax(const Kernel::Triangle_3& triangle) {
  return std::minmax({triangle.vertex(0).z(), triangle.vertex(1).z(), triangle.vertex(2).z()});
}

void help() {
  std::cout << "Usage: my_program [OPTIONS]\n"
            << "Options:\n"
            << "  --input <file>          Input file path\n"
            << "  --output <file>         Output file path\n"
            << "  --threshold <number>    Expansion threshold\n"
            << "  --help, -h              Show this help message\n";
}

std::string extract_base_filename(const std::string& path) {
  // Find the last directory separator
  size_t separator = path.find_last_of("/\\");
  std::string filename = (separator != std::string::npos)? path.substr(separator + 1) : path;
  return filename;
}

bool read_obj(
  const std::string input_file,
  std::vector<Kernel::Point_3> &vertices,
  std::vector<Face> &faces
) {
  // Read file
  std::ifstream input_stream;
  input_stream.open(input_file);
  if (!input_stream.is_open()) {
    std::cerr << "Failed to open input file: " << input_file << std::endl;
    return false;
  }
  std::string line;
  int current_material = -1;
  
  // Parse line by line
  while (getline(input_stream, line)) {
    try {
      std::istringstream line_stream(line);
      std::string line_type;
      line_stream >> line_type;
      
      // Vertex
      if (line_type == "v") {
        double x, y, z;
        line_stream >> x >> y >> z;
        vertices.emplace_back(Kernel::Point_3(x, y, z));
      }
      // Face
      if (line_type == "f") {
        faces.emplace_back();
        faces.back().id = faces.size()-1;
        faces.back().material = current_material;
        faces.back().block_id = -1;
        std::vector<Kernel::Point_3> face_vertices;
        int v;
        while (!line_stream.eof()) {
          line_stream >> v;
          face_vertices.emplace_back(vertices[v-1]);
        }
        if (face_vertices.size() == 3) {
          // introduces a triangle t with 3 vertices.
          faces.back().triangle = { face_vertices[0], face_vertices[1], face_vertices[2] };
          continue;
        }
      }
      // Material
      if (line_type == "usemtl") {
        line_stream >> current_material;
      }
    } catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
      return false;
    }
  }

  // Calculate and print max z
  double max_z = std::numeric_limits<double>::lowest();
  for (const auto& vertex : vertices) {
      if (vertex.z() > max_z) {
          max_z = vertex.z();
      }
  }

  return true;
}

void fan_triangulation_3(
  const std::vector<Kernel::Point_3> &convex_hull_3, 
  std::vector<std::array<int, 3>> &faces
) {
  int hull_size = convex_hull_3.size()/2;
  if (hull_size < 3) {
    std::cout << "Hull size < 3, not enough points for triangulation." << std::endl;
    return;
  }
  // top and bottom faces (horizontal)
  for (int i=1; i<hull_size-1; i++) {
    faces.push_back({0, i+1, i});
  }
  for (int i=hull_size+1; i<convex_hull_3.size()-1; i++) {
    faces.push_back({hull_size, i, i+1});
  }
  // wall faces (vertical)
  for (int i=0; i<hull_size-1; i++) {
    faces.push_back({i, i+1, hull_size+i});
    faces.push_back({i+1, hull_size+i+1, hull_size+i});
  }
  faces.push_back({hull_size-1, 0, hull_size*2-1});
  faces.push_back({0, hull_size, hull_size*2-1});
}

bool write_obj(
  const std::vector<Block> &blocks,
  const std::vector<std::vector<std::array<int, 3>>> &output_faces_chunk,
  const std::string output_file
) {
  std::ofstream out_stream;
  out_stream.open(output_file);

  if (!out_stream.is_open()) {
    std::cerr << "Failed to open output file: " << output_file << std::endl;
    return false;
  }
  // My name tag and mtl file
  out_stream << "# Created by MCHU\n" << "# Generated by MCHU\n";
  out_stream << "mtllib " << extract_base_filename(output_file) << ".mtl" << "\n";

  // Write vertices per block
  for (const auto& block: blocks) {
    for (const auto& vertex: block.convex_hull_3) {
      out_stream << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }
  }

  // Write faces per block
  int offset_amt = 1;
  for (int i = 0; i < blocks.size(); i++) {
    offset_amt += i==0 ? 0 : blocks[i-1].convex_hull_3.size();
    out_stream << "o Block " << i << "\n";
    for (const auto& triangle_idx : output_faces_chunk[i]) {
      out_stream << "usemtl " << i << "\n";
      out_stream << "f " 
        << triangle_idx[0] + offset_amt << " "
        << triangle_idx[1] + offset_amt << " "
        << triangle_idx[2] + offset_amt << "\n";
    }
  }
 
  out_stream.close();
  return true;
}

bool write_random_mtl(
  const int n,
  const std::string mtl_file
) {
  std::ofstream out_stream;
  out_stream.open(mtl_file);
  
  if (!out_stream.is_open()) {
    std::cerr << "Failed to open MTL file: " << mtl_file << std::endl;
    return false;
  }

  for (size_t i = 0; i < n; i++) {
    double r = std::round(((double) rand() / (RAND_MAX)) * 100.0f) / 100.0f;
    double g = std::round(((double) rand() / (RAND_MAX)) * 100.0f) / 100.0f;
    double b = std::round(((double) rand() / (RAND_MAX)) * 100.0f) / 100.0f;

    out_stream << "newmtl " << i << "\n";
    out_stream << "Ka 0.2 0.2 0.2\n"; // Ambient color
    out_stream << "Kd " << r << " " << g << " " << b << "\n"; // Diffuse color
    out_stream << "Ks 0.0 0.0 0.0\n"; // Specular color
    out_stream << "d 1.0\n";      // Transparency
    out_stream << "illum 2\n\n";  // Illumination model
  }

  out_stream.close();
  return true;
}

int main(int argc, const char *argv[]) {
  // Example: parse --input <file>, --output <file>, and --threshold <double>
  std::string input_file, output_file;
  double expansion_threshold = 1.0;
  if (argc==1) {
    std::cout << "Parse input and output .obj file path as parameters, for example:\n";
    help();
    return 0;
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      help();
      return 0;
    } else if ((arg == "--input" || arg == "-i") && i + 1 < argc) {
      input_file = argv[++i];
    } else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
      output_file = argv[++i];
    } else if ((arg == "--threshold" || arg == "-t") && i + 1 < argc) {
      expansion_threshold = atof(argv[++i]);
    } else {
      std::cout << "Parse input and output .obj file path as parameters, for example:\n";
      help();
      return 0;
    }
  }

  std::cout << "----------------------------------------------------------" << std::endl;
  std::cout << "Input file: " << input_file << std::endl;
  std::cout << "Output file: " << output_file << std::endl;
  std::cout << "Expansion threshold: " << expansion_threshold << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  // Main part begins here
  std::vector<Kernel::Point_3> input_vertices;
  std::vector<Face> input_faces;
  std::vector<Box> boxes;
  std::vector<Block> blocks;
  
  if (read_obj(input_file, input_vertices, input_faces)) {
    std::cout << "File loaded" << std::endl;
  } else {
    std::cout << "Input file loading failed" << std::endl;
    return -1;
  }
  
  // Generate expanded bounding boxes
  for (Faces_iterator current_face = input_faces.begin(); current_face != input_faces.end(); ++current_face) {
    CGAL::Bbox_3 face_box = current_face->triangle.bbox();
    CGAL::Bbox_2 expanded_box(
      face_box.xmin() - expansion_threshold, face_box.ymin() - expansion_threshold,
      face_box.xmax() + expansion_threshold, face_box.ymax() + expansion_threshold
    );
    boxes.push_back(Box(expanded_box, current_face));
  }
  
  // Compute neighbours
  std::vector<Faces_intersection> intersections;
  Box_intersector box_intersector(std::back_inserter(intersections));
  CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), box_intersector);
  std::cout << intersections.size() << " intersections found" << std::endl;
  for (auto const &intersection: intersections) {
    intersection.first->neighbours.emplace(intersection.second->id);
    intersection.second->neighbours.emplace(intersection.first->id);
  }
  
  // Assign blocks to faces: works like a simple growing algorithm
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
      }
      faces_in_current_block.pop_front();
    }
    ++current_block;
  }
  std::cout << current_block << " blocks found" << std::endl;
  
  // Create blocks and put points in each
  blocks.resize(current_block);
  for (auto &face: input_faces) {
    // Construct Point_2 right in the vector using emplace_back()
    blocks[face.block_id].points.emplace_back(face.triangle.vertex(0).x(), face.triangle.vertex(0).y());
    blocks[face.block_id].points.emplace_back(face.triangle.vertex(1).x(), face.triangle.vertex(1).y());
    blocks[face.block_id].points.emplace_back(face.triangle.vertex(2).x(), face.triangle.vertex(2).y());

    // Compute z value span (min/max)
    auto minmax = triangle_z_minmax(face.triangle);
    if (blocks[face.block_id].max < minmax.second) {
      blocks[face.block_id].max = minmax.second;
    }
    if (blocks[face.block_id].min > minmax.first) {
      blocks[face.block_id].min = minmax.first;
    }
  }
  
  // Compute convex hulls
  for (auto &block: blocks) {
    CGAL::convex_hull_2(block.points.begin(), block.points.end(), std::back_inserter(block.convex_hull));
    for (const auto& vertex: block.convex_hull) {
      block.convex_hull_3.emplace_back(vertex.x(), vertex.y(), block.min);
    }
    for (const auto& vertex: block.convex_hull) {
      block.convex_hull_3.emplace_back(vertex.x(), vertex.y(), block.max);
    }
  }

  // I think using indices would be faster
  // std::vector<Kernel::Triangle_3>> output_faces;
  std::vector<std::vector<std::array<int, 3>>> output_faces;
  output_faces.resize(current_block);
  for (int i = 0; i < current_block; i++) {
    fan_triangulation_3(blocks[i].convex_hull_3, output_faces[i]);
  }

  // Write output
  std::cout << "Writing files..." << std::endl;
  bool success = write_random_mtl(current_block, output_file+".mtl") && write_obj(blocks, output_faces, output_file);
  if (success) {
    std::cout << "Finished" << std::endl;
  } else {
    std::cout << "Output file writing failed" << std::endl;
  }
  return 0;
}
