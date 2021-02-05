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

	void CMesh::calJacobian(const std::size_t& ele, const std::size_t pg, Eigen::Matrix<double, 2, 2>& Jac, double& det)
	{
		Eigen::Matrix<double, 2, 4> ecoord;
		for (std::size_t i = 0; i < 4; i++)
		{
			auto nd = elements[ele].index_nd[i];
			ecoord(0,i) = nodes[nd].x;
			ecoord(1,i) = nodes[nd].y;
		}
		std::cout << ecoord << std::endl;
		Eigen::Vector2d lcoord;
		lcoord[0] = CQuadrature::qp_coords(pg,0);
		lcoord[1] = CQuadrature::qp_coords(pg, 1);
		Eigen::Matrix<double, 4, 2> spderiv = ShapeDeriv(lcoord);
		std::cout << lcoord << std::endl;
		std::cout << spderiv << std::endl;
		Jac = ecoord * spderiv;
		std::cout << "Jac=" << Jac << std::endl;
		auto inv = Jac.inverse();
		det = Jac.determinant();
	}
}

