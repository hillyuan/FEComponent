#include "control.hpp"

#include "yaml-cpp/yaml.h"
#include <string>
#include <filesystem>

namespace ROLLFEM2D
{
	CControl::CControl(char* file)
	{
		if (!std::filesystem::exists(file)) {
			std::cout << "Input file: " << file << " not found!\n" ;
			throw std::runtime_error("Mesh file not defined!");
		}

		YAML::Node config = YAML::LoadFile(file);
		if (YAML::Node doc = config["Mesh File"]) {
			const std::string filename = config["Mesh File"].as<std::string>();
			mesh.readin(filename.c_str());
		}
		else {
			throw std::runtime_error("Mesh file not defined!");
		}

		const YAML::Node& matls = config["Material"];
		for (YAML::const_iterator it = matls.begin(); it != matls.end(); ++it) {
			const YAML::Node& matls = *it;
			double youngs = matls["Youngs Modulus"].as<double>();
			double poisson = matls["Poisson Ratio"].as<double>();
			std::string mname = matls["Name"].as<std::string>();
			CMaterial matl(mname, youngs, poisson);
			mesh.materials.emplace_back(matl);
			matl.print(std::cout);
		}
	};

}

