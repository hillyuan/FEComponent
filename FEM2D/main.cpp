#include <fstream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <omp.h>

#include "mesh.hpp"
#include "control.hpp"


/**
 * \brief show available command line option
 */
void help() {
	printf("usage: rollfem2d [options] \n");
	printf(" -f: filename: control file name (neccesary)\n");
	printf(" -h: Show this help message. (optional)\n");
	printf(" -t <n>: Set number of OpenMP threads (optional)\n");
	exit(0);
}

void set_num_threads(char* arg) {
	int exec_threads;

	if (arg == NULL) {
		fprintf(stderr, "Error : specify number of OpenMP threads.\n");
		fprintf(stderr, "Format: -t <n>\n");
		exit(1);
	}

	exec_threads = atoi(arg);

	if (exec_threads == 0) {
		fprintf(stderr, "Error : specify 1 or more OpenMP threads.\n");
		exit(1);
	}
	omp_set_num_threads(exec_threads);
}



int main(int argc, char *argv[])
{
	char date[64];
	time_t t = time(NULL);
	std::cout << "=================================================================\n";
	std::cout << "=                                                               =\n";
	std::cout << "=                       RollFEM2D 1.0                           =\n";
	std::cout << "=                         March 2021                            =\n";
	std::cout << "=                                                               =\n";
	std::cout << "=================================================================\n\n";

	printf("Host info    \n");
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


	std::clock_t c_start = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();

	bool ok = false;
	ROLLFEM2D::CControl control;
	for (unsigned int i = 1; i < argc; i++) {
		if (strncmp("-f", argv[i], 2)==0) {
			control = ROLLFEM2D::CControl(argv[i+1]);
			ok = true; ++i;
		} else if (strncmp("-h", argv[i], 2) == 0) {
			help();
		}
		else if (strncmp("-t", argv[i], 2) == 0) {
			set_num_threads(argv[i+1]);
		}
	}
	if (!ok) {
		help();
		return -1;
	}
	printf("  threads:    %d\n\n", omp_get_max_threads());

	//ROLLFEM2D::CControl control(argv[1]);

	printf("Calculating......  \n");
	control.calGlobalStiffMatrix();
	control.ApplyDistributedLoads();
	control.ApplyConstraints();
	control.Solve();
	std::clock_t c_end = std::clock();
	std::cout << "  calculate time:  " << (c_end - c_start) / CLOCKS_PER_SEC << " s\n\n";

	printf("Outputing......  \n\n");
	control.Update();
	int c = control.VTKOutput();
	c = control.CSVOutput();
	
	c_end = std::clock();
	auto t_end = std::chrono::high_resolution_clock::now();

	printf("Summary:  \n");
	std::cout << std::fixed << std::setprecision(2) << "  CPU time used: "
		<< (c_end - c_start) / CLOCKS_PER_SEC << " s\n"
		<< "  Wall clock time passed: "
		<< std::chrono::duration<double>(t_end - t_start).count()
		<< " s\n";
}
