#include <fstream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "mesh.hpp"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "Usage: rollfem2d filename(vtk file)";
		return -1; 
	}
	std::clock_t c_start = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();

	ROLLFEM2D::CMesh mesh;
	mesh.readin(argv[1]);
	
	std::clock_t c_end = std::clock();
	auto t_end = std::chrono::high_resolution_clock::now();

	std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
		<< (c_end - c_start) / CLOCKS_PER_SEC << " s\n"
		<< "Wall clock time passed: "
		<< std::chrono::duration<double>(t_end - t_start).count()
		<< " s\n";
}