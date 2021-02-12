#include "mesh.hpp"

#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	int CMesh::readin(char* meshfile)
	{
		std::ifstream input(meshfile, std::ios_base::in);
		if (input.fail()) {
			std::cout << "Mesh file not found!\n";
			return -2;
		}

		std::string line;
		std::string dummy, name, dtype;
		std::istringstream tokenizer;

		double dummyz;

		while (!input.eof()) {
			if (line.find("DATASET") == 0) {
				if (line.find("UNSTRUCTURED_GRID") == std::string::npos) {
					std::cout << "Mesh format not supported!\n";
					input.close();
					return -2;
				}
				std::getline(input, line);
			}
			else if (line.find("POINTS") == 0) {
				std::cout << line << std::endl;
				tokenizer.str(line.substr(0, line.find_last_not_of(" \r\n") + 1));
				tokenizer.clear();
				tokenizer >> dummy >> num_nodes;
				this->nodes.resize(num_nodes);
				for (int i = 0; i < num_nodes; ++i) {
					input >> nodes[i].x >> nodes[i].y >> dummyz;
				}
				std::getline(input, line);
			}
			else if (line.find("CELLS") == 0) {
				std::cout << line << std::endl;
				tokenizer.str(line.substr(0, line.find_last_not_of(" \r\n") + 1));
				tokenizer.clear();
				tokenizer >> dummy >> num_elements >> dummy;
				this->elements.resize(num_elements);
				for (size_t i = 0; i < num_elements; ++i) {
					std::getline(input, line);
					tokenizer.str(line.substr(0, line.find_last_not_of(" \r\n") + 1));
					tokenizer.clear();
					tokenizer >> dummy >> elements[i].index_nd[0] >> elements[i].index_nd[1]
					          >> elements[i].index_nd[2] >> elements[i].index_nd[3];
				}
				std::getline(input, line);
			}
			else {
				std::getline(input, line);
			}
		}
		input.close();

		//print_elements(std::cout);

		return 0;
	}

	void CMesh::print_elements(std::ostream& os) const
	{
		os << "Number of elements = " << num_elements << std::endl;
		os << "Element ID:    Nodes List\n";
		std::size_t ne = -1;
		for ( auto ele: elements )
		{
			os << ++ne << "  " << ele.index_nd[0] << "  " << ele.index_nd[1]
			   << "  " << ele.index_nd[2] << "  " << ele.index_nd[3] << std::endl;
		}
	}

	void CMesh::print_nodes(std::ostream& os) const
	{
		os << "Number of nodes = " << num_nodes << std::endl;
		os << "Node ID:    Coordinate\n";
		std::size_t cnt = -1;
		for ( auto nd: nodes)
		{
			os << ++cnt << "  " << nd.x  <<  "  " << nd.y << std::endl;
		}
	}

	void CMesh::calElementalStiffMatrix(const std::size_t& ele, Eigen::Matrix<double, 8, 8>& K
		, std::vector<T>& triplets)
	{
		Eigen::Matrix<double, 2, 4> ecoord;
		Eigen::Matrix<double, 2, 2> Jac;
		Eigen::Matrix<double, 3, 8> B;
		int imatl = elements[ele].matl_id;

		std::size_t n=3;
		for (std::size_t i = 0; i<4; i++ )
		{
			auto nd = elements[ele].index_nd[i];
			ecoord(0,n) = nodes[nd].x;
			ecoord(1,n) = nodes[nd].y;
			n--;
		}
	//	std::cout << "ecood:" << ecoord << std::endl;

		Eigen::Vector2d lcoord;
		K = Eigen::Matrix<double, 8, 8>::Zero();
		for (std::size_t npg = 0; npg < 4; npg++)
		{
			Eigen::Matrix<double, 4, 2> spderiv = CQuadrature::ShapeDerivs[npg];
			Jac = ecoord * spderiv;
			auto inv = Jac.inverse();
			double wg = CQuadrature::weights[npg] * Jac.determinant();
			auto gderiv = spderiv * inv;
			for (std::size_t i = 0; i < 4; i++)
			{
				B(0, 2 * i) = gderiv(i, 0);
				B(0, 2 * i + 1) = 0.0;
				B(1, 2 * i) = 0.0;
				B(1, 2 * i + 1) = gderiv(i, 1);
				B(2, 2 * i) = gderiv(i, 1);
				B(2, 2 * i + 1) = gderiv(i, 0);
			}
			K += wg * B.transpose() * materials[imatl].ElasticMatrix * B;
		}

	//	std::cout << K << std::endl;
		for (std::size_t i = 0; i < 4; i++)
		{
			for (std::size_t j = 0; j < 4; j++)
			{
				T trplt11(2 * elements[ele].index_nd[i] + 0, 2 * elements[ele].index_nd[j] + 0, K(2 * i + 0, 2 * j + 0));
				T trplt12(2 * elements[ele].index_nd[i] + 0, 2 * elements[ele].index_nd[j] + 1, K(2 * i + 0, 2 * j + 1));
				T trplt21(2 * elements[ele].index_nd[i] + 1, 2 * elements[ele].index_nd[j] + 0, K(2 * i + 1, 2 * j + 0));
				T trplt22(2 * elements[ele].index_nd[i] + 1, 2 * elements[ele].index_nd[j] + 1, K(2 * i + 1, 2 * j + 1));

				triplets.push_back(trplt11);
				triplets.push_back(trplt12);
				triplets.push_back(trplt21);
				triplets.push_back(trplt22);
			}
		}
	}

	void CMesh::calGlobalStiffMatrix()
	{
		std::vector<T> tripletList;
		tripletList.reserve(num_elements*2*10);

		std::size_t num_dofs = 2 * num_elements;
		StiffMatrix = Eigen::SparseMatrix<double>(num_dofs, num_dofs);

		Eigen::Matrix<double, 8, 8> eleK;
		for( std::size_t i=0; i< num_elements; i++ )
		{ 
			this->calElementalStiffMatrix(i, eleK, tripletList);
		}
		StiffMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	}
}

