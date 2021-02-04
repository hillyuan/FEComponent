#ifndef ROLLFEM2D_MESH_HPP
#define ROLLFEM2D_MESH_HPP

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <Eigen\Eigen>

namespace ROLLFEM2D
{
	struct CNode
	{
		double x, y;
	};

	struct CElement
	{
		std::size_t n0, n1, n2, n3;
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
	};
}


#endif
