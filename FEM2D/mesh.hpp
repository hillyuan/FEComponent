#ifndef ROLLFEM2D_MESH_HPP
#define ROLLFEM2D_MESH_HPP

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <Eigen\Eigen>

#include "material.hpp"
#include "element.hpp"

namespace ROLLFEM2D
{
	struct CNode
	{
		double x, y;
	};

	struct CEdge
	{
		std::size_t id_element;   // index of element
		std::size_t n_edge;       // edge number of the element

		CEdge(std::size_t& ie, std::size_t& nd)
			: id_element(ie), n_edge(nd)
		{};
	};

	class CMesh
	{
		typedef Eigen::Triplet<double> T;

	public:
		std::size_t num_nodes;
		std::size_t num_elements;
		std::vector<CNode>         nodes;
		std::vector<CElement>      elements;
		std::vector<CMaterial>     materials;

		std::map< std::string, std::vector<std::size_t> > NodeSets;
		std::map< std::string, std::vector<std::size_t> > ElementSets;
		std::map< std::string, std::vector<CEdge> >       SideSets;

		int readin(const char *);
		void print_elements(std::ostream& os) const;
		void print_nodes(std::ostream& os) const;
		Eigen::Vector2d getCenterCoord(const std::size_t& ele) const;

		void calElementalStiffMatrix (const std::size_t&, std::array<T,64> & triplets);
		void calGlobalStiffMatrix(Eigen::SparseMatrix<double>&);
		void updateElements(Eigen::VectorXd&);

	private:
		std::size_t n_wele, n_iele, n_bele;
		std::vector<double> ThicknessDefine;

		void ConstructElementThickness();
	};
}


#endif
