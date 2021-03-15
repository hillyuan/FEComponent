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

		int pbtype = 0;
		if (YAML::Node mtype = config["Problem Type"]) {
			std::string tt = mtype.as<std::string>();
			if (tt == "Plane Strain") pbtype = 1;
		}

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
				CMaterial matl(pbtype, mname, youngs, poisson);
				mesh.materials.emplace_back(matl);
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
		if (YAML::Node doc = config["VTK File"]) {
			outfile = doc.as<std::string>();
		}

		if (YAML::Node csvdef = config["CSV File"]) {
			for (YAML::const_iterator it = csvdef.begin(); it != csvdef.end(); ++it) {
				const YAML::Node& dl = *it;
				auto fname = dl["File Name"].as<std::string>();
				auto outype = dl["TYPE"].as<std::string>();
				auto nset = dl["NSET"].as<std::string>();
				if (mesh.NodeSets.find(nset) == mesh.NodeSets.end())
				{
					std::cout << "NodeSet:" << nset << " not found. CSV File ignored!\n";
				}
				else {
					cvsout myout(fname, nset, outype);
					cvsouts.emplace_back(myout);
				}
			}
		}

		loads.resize(2 * mesh.num_nodes);
		loads.setZero();
		forces.resize(2 * mesh.num_nodes);
		forces.setZero();
	};


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
			//	auto force = stress.transpose()*ele.B;
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

	void CControl::calEquivalentNodalForce()
	{
		Eigen::Vector<double, 8> eforce;
//#pragma omp parallel for private(force)
		for (int i = 0; i < mesh.num_elements; ++i)
		{
			eforce = mesh.elements[i].calNodalForce();
			for (int j=0;j<4;++j)
			{
				std::size_t nd = mesh.elements[i].index_nd[j];
//#pragma omp critical
				forces(2 * nd) += eforce(j * 2);
				forces(2 * nd + 1) += eforce(j * 2 + 1);
			}
		}
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
		output << "VECTORS Force double" << std::endl;
		cnt = -1;
		for (int i = 0; i < mesh.num_nodes; i++) {
			output << forces[++cnt] << " " << forces[++cnt] << " 0.0" << std::endl;
		}

		output << "CELL_DATA " << mesh.num_elements << std::endl;
		output << "SCALARS Strain double 3" << std::endl;
		output << "LOOKUP_TABLE default" << std::endl;
		Eigen::Vector3d ss;
		for (int i = 0; i < mesh.num_elements; ++i) {
			ss.setZero();
			for (int j = 0; j < 4; ++j) {
				ss += mesh.elements[i].strain[j];
			}
			ss *= 0.25;
			output << ss[0] << " " << ss[1] << " " << ss[2] << std::endl;
		}
		output << "SCALARS Stress double 3" << std::endl;
		output << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i < mesh.num_elements; ++i) {
			ss.setZero();
			for (int j = 0; j < 4; ++j) {
				ss += mesh.elements[i].stress[j];
			}
			ss *= 0.25;
			output << ss[0] << " " << ss[1] << " " << ss[2] << std::endl;
		}

		output.close();
		return 0;
	}

	int CControl::CSVOutput() const
	{
		for (auto csv : cvsouts) {
			std::ofstream output(csv.fname.c_str(), std::ios_base::out);
			if (output.fail()) {
				return -2;
			}
			std::vector<std::size_t> outnodes = mesh.NodeSets.at(csv.ndset);
			output << "#Output of " << csv.type << " of nodeset " << csv.ndset << std::endl;
			if (csv.type == "y-disp") {
				for (auto nd : outnodes) {
					output << mesh.nodes[nd].x << "," << displacements(nd * 2 + 1) << std::endl;
				}
			} else if (csv.type == "y-force") {
				for (auto nd : outnodes) {
					output << mesh.nodes[nd].x << "," << forces(nd * 2 + 1) << std::endl;
				}
			}
			output.close();
		}

		return 0;
	}

}

