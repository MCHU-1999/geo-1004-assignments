cmake_minimum_required(VERSION 3.12)

project(merge_buildings)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core)

# Create executable
add_executable(${PROJECT_NAME} hw1.cpp)

# Link CGAL to your executable
target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::CGAL_Core)

#end of the file