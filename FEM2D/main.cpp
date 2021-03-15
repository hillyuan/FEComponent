#include <fstream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <omp.h>


#include "mesh.hpp"
#include "control.hpp"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cout << "Usage: rollfem2d filename(control file)";
		return -1; 
	}

	char date[64];
	time_t t = time(NULL);
	printf("execute......  \n");
	char* libvar;
#ifdef _WINDOWS
	libvar = getenv("COMPUTERNAME");
	printf("  host:       %s\n", libvar);
	libvar = getenv("USERNAME");
	printf("  user:       %s\n", libvar);
#else
	libvar = getenv("HOSTNAME");
	printf("  host:       %s\n", libvar);
	libvar = getenv("USER");
	printf("  user:       %s\n", libvar);
#endif
	strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S%z", localtime(&t));
	printf("  date:       %s\n", date);
	printf("  threads:    %d\n\n", omp_get_max_threads());


	std::clock_t c_start = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();

	ROLLFEM2D::CControl control(argv[1]);

	printf("calculating......  \n");
	control.calGlobalStiffMatrix();
	control.ApplyDistributedLoads();
	control.ApplyConstraints();
	control.Solve();
	std::clock_t c_end = std::clock();
	std::cout << "  calculate time:  " << (c_end - c_start) / CLOCKS_PER_SEC << " sec\n\n";

	printf("outputing......  \n\n");
	control.Update();
	int c = control.VTKOutput();
	c = control.CSVOutput();
	
	c_end = std::clock();
	auto t_end = std::chrono::high_resolution_clock::now();

	printf("summary:  \n");
	std::cout << std::fixed << std::setprecision(2) << "  CPU time used: "
		<< (c_end - c_start) / CLOCKS_PER_SEC << " s\n"
		<< "  Wall clock time passed: "
		<< std::chrono::duration<double>(t_end - t_start).count()
		<< " s\n";
}