cmake_minimum_required(VERSION 3.14)
project(ROLLFEM2D CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(YamlCpp REQUIRED)
if(YAMLCPP_FOUND)
  message(STATUS "Yaml-cpp Found. ${YAMLCPP_INCLUDE_DIR}")
  include_directories(${YAMLCPP_INCLUDE_DIR})
  link_directories(${YAMLCPP_LIBRARY_DIR})
endif()

find_package(Eigen REQUIRED)
if(EIGEN_FOUND)
  message(STATUS "Eigen ${EIGEN_VERSION} Found. ${EIGEN_INCLUDE_DIR}")
  include_directories(${EIGEN_INCLUDE_DIR})
endif()

file(GLOB SOURCE src/*.cpp src/*/*.cpp)
file(GLOB HEADER src/*.hpp src/*/*.hpp
     ${CMAKE_CURRENT_BINARY_DIR}/include/*.hpp)

add_executable(rollfem2d main.cpp control.cpp mesh.cpp material.cpp element.cpp)
target_link_libraries(rollfem2d debug ${YAMLCPP_LIBRARY_DEBUG})
target_link_libraries(rollfem2d optimized ${YAMLCPP_LIBRARY_RELEASE})
