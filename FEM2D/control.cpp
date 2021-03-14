#include "control.hpp"

#include "yaml-cpp/yaml.h"
#include <string>
#include <filesystem>
#include <cassert>

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

		char* p;
		if (YAML::Node bnd = config["Constraint"]) {
			for (YAML::const_iterator it = bnd.begin(); it != bnd.end(); ++it) {
				const YAML::Node& bnd = *it;
				int uxy = bnd["Type"].as<int>();
				double val = bnd["Value"].as<double>();
				std::string mname = bnd["NSET"].as<std::string>();
				int num = strtol(mname.c_str(), &p, 10);
				if(*p == 0) { 
					Constraint cst(num, uxy, val);
					constraints.emplace_back(cst);
				}
				else {
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
		}
		else {
			throw std::runtime_error("Constraint condition not defined!");
		}

		if (YAML::Node dload = config["DLoad"]) {
			for (YAML::const_iterator it = dload.begin(); it != dload.end(); ++it) {
				const YAML::Node& dl = *it;
				std::vector<double> val = dl["Value"].as<std::vector<double>>();
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

		outfile = "out.vtk";
		if (YAML::Node doc = config["Output File"]) {
			outfile = doc.as<std::string>();
		}
		if (YAML::Node doc = config["Thickness Define"]) {
			ThicknessDefine = doc.as<std::vector<double>>();
		}

		loads.resize(2 * mesh.num_nodes);
		loads.setZero();
	};

	void CControl::ConstructElementThickness()
	{
		assert(ThicknessDefine.size() == 11);
		double sy = ThicknessDefine[0];
		double wD1 = 0.5*ThicknessDefine[1];
		double wL1 = 0.5*ThicknessDefine[2];
		double wD2 = 0.5*ThicknessDefine[3];
		double iD1 = 0.5*ThicknessDefine[4];
		double iL1 = 0.5*ThicknessDefine[5];
		double iD2 = 0.5*ThicknessDefine[6];
		double offset = ThicknessDefine[7];
		double bD1 = 0.5*ThicknessDefine[8];
		double bL1 = 0.5*ThicknessDefine[9];
		double bD2 = 0.5*ThicknessDefine[10];
		// inter roll
		double dy;
		Eigen::Vector2d center;
		double cy = sy + wD1 + iD1;
		for (int i = 0; i < mesh.n_iele; ++i)
		{
			center = mesh.getCenterCoord(i);
			dy = center[1] - cy;
			if( center[0]>-iL1+ offset && center[0] < iL1+ offset) // with diameter D1
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(iD1 * iD1 - dy * dy);
			} else // with diameter D3
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(iD2 * iD2 - dy * dy);
			}
		}

		// back roll
		cy = cy + iD1 + bD1;
		for (int i = mesh.n_iele; i < mesh.n_bele; ++i)
		{
			center = mesh.getCenterCoord(i);
			dy = center[1] - cy;
			if (center[0] > -bL1 && center[0] < bL1) // with diameter D1
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(bD1 * bD1 - dy * dy);
			}
			else // with diameter D3
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(bD2 * bD2 - dy * dy);
			}
		}

		// work roll
		for (int i = mesh.n_bele; i < mesh.n_wele; ++i)
		{
			center = mesh.getCenterCoord(i);
			dy = center[1] - cy;
			if (center[0] > -wL1 && center[0] < wL1) // with diameter D1
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(wD1 * wD1 - dy * dy);
			}
			else // with diameter D3
			{
				mesh.elements[i].thick = 2.0 * std::sqrt(wD2 * wD2 - dy * dy);
			}
		}

	}

	void CControl::ApplyConstraints()
	{
		std::vector<std::size_t> indicesToConstraint;

		for (auto cst : constraints)
		{
			if (cst.id_node >= 0) {
				if (cst.type & Constraint::UX)
				{
					indicesToConstraint.push_back(2 * cst.id_node);
				}
				if (cst.type & Constraint::UY)
				{
					indicesToConstraint.push_back(2 * cst.id_node + 1);
				}
			}
			else {
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

		for (auto idit : indicesToConstraint)
		{
			loads(idit) = 0.0;
		}

	//	std::cout << StiffMatrix << std::endl;
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
				normal[0] = 0.5 * (mesh.nodes[nd1].y - mesh.nodes[nd0].y);
				normal[1] = 0.5 * (mesh.nodes[nd0].x - mesh.nodes[nd1].x);
				if (load.type == 2) {   // surface presuure
					double area = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
					normal[0] = area * load.direction[0];
					normal[1] = area * load.direction[1];
				}
				loads(2 * nd0) += load.val * normal[0];
				loads(2 * nd0 + 1) += load.val * normal[1];
				loads(2 * nd1) += load.val * normal[0];
				loads(2 * nd1 + 1) += load.val * normal[1];
				
			}
		}

//		std::cout << loads << std::endl;
	}

	void CControl::ApplyInitialStrain()
	{
		Eigen::Matrix<double, 3, 3> D = mesh.materials[0].ElasticMatrix;
		for (auto load : initstrain)
		{
			// consider edge pressure only
			auto eleids = mesh.ElementSets[load.SetName];
			for (int i = 0; i < eleids.size(); ++i)
			{
				CElement ele = mesh.elements[eleids[i]];
				Eigen::Vector3d stress = D * load.strain;
				auto force = stress.transpose()*ele.B;
			//	loads(2 * nd0) += load.val * normal[0];
			//	loads(2 * nd0 + 1) += load.val * normal[1];
			//	loads(2 * nd1) += load.val * normal[0];
			//	loads(2 * nd1 + 1) += load.val * normal[1];

			}
		}

		//		std::cout << loads << std::endl;
	}

	void CControl::Solve()
	{
		Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
		solver.compute(StiffMatrix);
		displacements = solver.solve(loads);

	//	std::cout << "displacament:=" << displacements << std::endl;
	}

	int CControl::VTKOutput() const
	{
		std::ofstream output(outfile.c_str(), std::ios_base::out);
		if (output.fail()) {
			return -2;
		}
		output << "# vtk DataFile Version 4.0" << std::endl;
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

		output << "CELL_DATA " << mesh.num_elements << std::endl;
		output << "SCALARS Strain double 3" << std::endl;
		output << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i < mesh.num_elements; ++i) {
			output << mesh.elements[i].strain[0][0] << " " << mesh.elements[i].strain[0][1] 
				<< " " << mesh.elements[i].strain[0][2] << std::endl;
		}
		output << "SCALARS Stress double 3" << std::endl;
		output << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i < mesh.num_elements; ++i) {
			output << mesh.elements[i].stress[0][0] << " " << mesh.elements[i].stress[0][1]
				<< " " << mesh.elements[i].stress[0][2] << std::endl;
		}

		output.close();
		return 0;
	}

}

