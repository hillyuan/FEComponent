#include "control.hpp"

#include "yaml-cpp/yaml.h"
#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	CControl::CControl(std::string& file)
	{
		YAML::Node config = YAML::LoadFile(file);
		//	std::ifstream fin(contolfile);
		//	if( !fin.good() )
		//		throw std::runtime_error("File not exist");
		//	YAML::Parser parser(fin);
		//	parser.HandleNextDocument(config);
	};

}

