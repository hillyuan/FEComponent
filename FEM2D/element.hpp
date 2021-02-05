#ifndef _ROLLFEM2D_ELEMENT_HPP
#define _ROLLFEM2D_ELEMENT_HPP

#include <iostream>
#include <valarray>
#include <string>
#include <Eigen\Eigen>

namespace ROLLFEM2D
{
	void ShapeFunc(double lcoord[2], double sfunc[4]);
	Eigen::Matrix<double, 4, 2> ShapeDeriv(Eigen::Vector2d& lcoord);

	Eigen::Matrix<double, 4, 2> make_quadrature();
	Eigen::Vector4d make_weights();

	struct CQuadrature
	{
		const static Eigen::Vector4d weights;
		const static Eigen::Matrix<double, 4, 2> qp_coords;
	};

	class CElement
	{
	public:
		int matl_id;
		std::size_t index_nd[4];
	};
}


#endif
