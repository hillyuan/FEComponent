#include "control.hpp"

#include "yaml-cpp/yaml.h"
#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	CControl::CControl(char* file)
	{
		YAML::Node config = YAML::LoadFile(file);
		const std::string username = config["Mesh File"].as<std::string>();
		const YAML::Node& matls = config["Material"];
		for (YAML::const_iterator it = matls.begin(); it != matls.end(); ++it) {
			const YAML::Node& matls = *it;
			double youngs = matls["Youngs Modulus"].as<double>();
			double poisson = matls["Poisson Ratio"].as<double>();
			std::string mname = matls["Name"].as<std::string>();
			CMaterial matl(mname, youngs, poisson);
			materials.emplace_back(matl);
			matl.print(std::cout);
		}
	};

}

