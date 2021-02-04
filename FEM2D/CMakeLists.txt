cmake_minimum_required(VERSION 3.14)
project(ROLLFEM2D CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# For OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Yaml-cpp REQUIRED)
if(Yaml-cpp_FOUND)
  message(STATUS "Yaml-cpp Found. ${YAML_CPP_INCLUDE_PATH}")
  include_directories(${YAML_CPP_INCLUDE_PATH})
endif()


file(GLOB SOURCE src/*.cpp src/*/*.cpp)
file(GLOB HEADER src/*.hpp src/*/*.hpp
     ${CMAKE_CURRENT_BINARY_DIR}/include/*.hpp)

add_executable(rollfem2d main.cpp)

