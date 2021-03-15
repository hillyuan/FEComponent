#ifndef _ROLLFEM2D_ELEMENT_HPP
#define _ROLLFEM2D_ELEMENT_HPP

#include <iostream>
#include <array>
#include <string>
#include <Eigen\Eigen>

namespace ROLLFEM2D
{
	Eigen::Vector2d ShapeFunc(Eigen::Vector2d& lcoord);
	Eigen::Matrix<double, 4, 2> ShapeDeriv(Eigen::Vector2d& lcoord);

	Eigen::Matrix<double, 4, 2> make_quadrature();
	Eigen::Vector4d make_weights();
	std::array<Eigen::Matrix<double, 4, 2>, 4> make_ShapeDerivs();

	struct CQuadrature
	{
		const static Eigen::Vector4d weights;
		const static Eigen::Matrix<double, 4, 2> qp_coords;
		const static std::array<Eigen::Matrix<double, 4, 2>,4> ShapeDerivs;
	};

	class CElement
	{
	public:
		int matl_id;
		double thick;
		std::size_t index_nd[4];
		double wg[4];

		Eigen::Matrix<double, 3, 8> B[4];

		Eigen::Vector3d strain[4];
		Eigen::Vector3d stress[4];

		Eigen::Vector<double, 8> calNodalForce() const;
	};
}


#endif
