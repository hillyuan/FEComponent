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

		if (YAML::Node matls = config["Material"]) {
			for (YAML::const_iterator it = matls.begin(); it != matls.end(); ++it) {
				const YAML::Node& matls = *it;
				double youngs = matls["Youngs Modulus"].as<double>();
				double poisson = matls["Poisson Ratio"].as<double>();
				std::string mname = matls["Name"].as<std::string>();
				CMaterial matl(mname, youngs, poisson);
				mesh.materials.emplace_back(matl);
				matl.print(std::cout);
			}
		}
		else {
			throw std::runtime_error("Material properties not defined!");
		}

		if (YAML::Node bnd = config["Constraint"]) {
			for (YAML::const_iterator it = bnd.begin(); it != bnd.end(); ++it) {
				const YAML::Node& bnd = *it;
				int uxy = bnd["Type"].as<int>();
				double val = bnd["Value"].as<double>();
				std::string mname = bnd["NSET"].as<std::string>();
				if (mesh.NodeSets.find(mname) == mesh.NodeSets.end())
				{
					std::cout << "NodeSet:" << mname << " not found. Constraint ignored!\n";
				}
				else {
					Constraint cst(mname, uxy, val);
					constraints.emplace_back(cst);
				}
			}
		}
		else {
			throw std::runtime_error("Constraint condition not defined!");
		}

		loads.resize(2 * mesh.num_nodes);
		loads.setZero();
	};

	void CControl::ApplyConstraints()
	{
		std::vector<std::size_t> indicesToConstraint;

		for (auto cst : constraints)
		{
			auto nodes = mesh.NodeSets[cst.NSetName];
			if (cst.type & Constraint::UX)
			{
				for (int j = 0; j < nodes.size(); ++j)
					indicesToConstraint.push_back(2 * nodes[j]);
			}
			if (cst.type & Constraint::UY)
			{
				for (int j = 0; j < nodes.size(); ++j)
					indicesToConstraint.push_back(2 * nodes[j] + 1);
			}
		}

		for (int k = 0; k < StiffMatrix.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(StiffMatrix, k); it; ++it)
			{
				for (auto idit : indicesToConstraint)
				{
					if (it.row() == idit || it.col() == idit)
					{
						it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
					}
				}
			}
		}

		std::cout << StiffMatrix << std::endl;
	}

}

