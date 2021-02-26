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

		if (YAML::Node dload = config["DLoad"]) {
			for (YAML::const_iterator it = dload.begin(); it != dload.end(); ++it) {
				const YAML::Node& dl = *it;
				double val = dl["Value"].as<double>();
				std::string mname = dl["SSET"].as<std::string>();
				if (mesh.SideSets.find(mname) == mesh.SideSets.end())
				{
					std::cout << "SideSet:" << mname << " not found. DLoad ignored!\n";
				}
				else {
					DLoad cst(mname, val);
					dloads.emplace_back(cst);
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

	void CControl::ApplyDistributedLoads()
	{
		std::size_t nd0, nd1, ie;
		double normal[2];
		for( auto load: dloads)
		{
			// consider edge pressure only
			auto edges = mesh.SideSets[load.SetName];
			for (int i = 0; i < edges.size(); ++i)
			{
				nd0 = edges[i].n_edge;
				nd1 = edges[i].n_edge == 3 ? 0 : edges[i].n_edge + 1;
				ie = edges[i].id_element;
				nd0 = mesh.elements[ie].index_nd[nd0];
				nd1 = mesh.elements[ie].index_nd[nd1];
				normal[0] = mesh.nodes[nd0].y - mesh.nodes[nd1].y;
				normal[1] = mesh.nodes[nd1].x - mesh.nodes[nd0].x;
				loads(2 * nd0) += load.val * normal[0];
				loads(2 * nd0 + 1) += load.val * normal[1];
				loads(2 * nd1) += load.val * normal[0];
				loads(2 * nd1 + 1) += load.val * normal[1];
			}
		}

		std::cout << loads << std::endl;
	}

	void CControl::solve()
	{
		Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
		solver.compute(StiffMatrix);
		displacements = solver.solve(loads);

		std::cout << displacements << std::endl;
	}

	int CControl::VTKOutput(std::string filename) const
	{
		std::ofstream output(filename.c_str(), std::ios_base::out);
		if (output.fail()) {
			return -2;
		}
		output << "# vtk DataFile Version 3.0" << std::endl;
		output << "Result of 3-Roller FEM" << std::endl;
		output << "ASCII" << std::endl;
		output << "DATASET UNSTRUCTURED_GRID" << std::endl;

		output << "POINTS " << mesh.num_nodes << " double" << std::endl;
		for (auto p: mesh.nodes) {
			output << std::setprecision(12) << p.x << " " << p.y << " 0.0" << std::endl;
		}

		output << "CELLS " << mesh.num_elements << " " << mesh.num_elements * 5 << std::endl;
		for (auto ele: mesh.elements) {
			output << 4 << " " << ele.index_nd[0] << " " << ele.index_nd[1] <<
				" " << ele.index_nd[2] << " " << ele.index_nd[3] << std::endl;
		}
		output << "CELL_TYPES " << mesh.num_elements << std::endl;
		for (int i = 0; i < mesh.num_elements; ++i) {
			output << 9 << std::endl;
		}

		output << "POINT_DATA " << mesh.num_nodes << std::endl;
		output << "VECTORS Displacement double" << std::endl;
		std::size_t cnt = -1;
		for (int i = 0; i < mesh.num_nodes; i++) {
			output << displacements[++cnt] << " " << displacements[++cnt] << " 0.0" << std::endl;
		}

		output.close();
		return 0;
	}

}

