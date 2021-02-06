#ifndef ROLLFEM2D_MESH_HPP
#define ROLLFEM2D_MESH_HPP

#include <iostream>
#include <vector>
#include <string>
#include <Eigen\Eigen>

#include "material.hpp"
#include "element.hpp"

namespace ROLLFEM2D
{
	struct CNode
	{
		double x, y;
	};

	struct CMesh
	{
		std::size_t num_nodes;
		std::size_t num_elements;
		std::vector<CNode>         nodes;
		std::vector<CElement>      elements;
		std::vector<CMaterial>     materials;

		Eigen::SparseMatrix<double> StiffMatrix;

		int readin(char *);
		void print_elements(std::ostream& os) const;
		void print_nodes(std::ostream& os) const;

		void calElementalStiffMatrix (const std::size_t&, Eigen::Matrix<double, 8, 8>&);
		void calGlobalStiffMatrix();
	};
}


#endif
