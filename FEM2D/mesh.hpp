#ifndef ROLLFEM2D_MESH_HPP
#define ROLLFEM2D_MESH_HPP

#include <iostream>
#include <valarray>
#include <string>

namespace ROLLFEM2D
{
	struct CMesh
	{
		std::size_t num_nodes;
		std::size_t num_elements;
		std::valarray<double>      coords;
		std::valarray<std::size_t> elements;

		int readin(char *);
		void print_elements(std::ostream& os) const;
		void print_nodes(std::ostream& os) const;
	};
}


#endif
