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

		int readin(char *);
		void print_elements(std::ostream& os) const;
		void print_nodes(std::ostream& os) const;

		void calJacobian(std::size_t&, std::size_t, Eigen::Matrix<double,4,2>&, double&) const;
	};
}


#endif
