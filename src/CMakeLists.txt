cmake_minimum_required(VERSION 2.8)
find_package(Eigen3 REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${EIGEN3_INCLUDE_DIR}")
add_definitions("-Wall -std=c++11") 
add_library(a STATIC
  Material.cpp
)
