#include "element.hpp"


namespace ROLLFEM2D
{
	void ShapeFunc(double lcoord[2], double sfunc[4])
	{
		sfunc[0] = 0.25 * (1.0 - lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[1] = 0.25 * (1.0 + lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[2] = 0.25 * (1.0 + lcoord[0]) * (1.0 + lcoord[1]);
		sfunc[3] = 0.25 * (1.0 - lcoord[0]) * (1.0 + lcoord[1]);
	}

	void ShapeDeriv(double lcoord[2], Eigen::Matrix<double, 3, 3>& spderiv)
	{
		spderiv(0,0) = -0.25 * (1.0 - lcoord[1]);
		spderiv(1,0) = 0.25 * (1.0 - lcoord[1]);
		spderiv(2,0) = 0.25 * (1.0 + lcoord[1]);
		spderiv(3,0) = -0.25 * (1.0 + lcoord[1]);

		spderiv(0,1) = -0.25 * (1.0 - lcoord[0]);
		spderiv(1,1) = 0.25 * (1.0 - lcoord[0]);
		spderiv(2,1) = 0.25 * (1.0 + lcoord[0]);
		spderiv(3,1) = -0.25 * (1.0 + lcoord[0]);
	}

	Eigen::Matrix<double, 4, 2> make_quadrature()
	{
		Eigen::Matrix<double, 4, 2> qp_coords;
		qp_coords(0, 0) = -0.577350269189626;
		qp_coords(0, 1) = -0.577350269189626;
		qp_coords(1, 0) = 0.577350269189626;
		qp_coords(1, 1) = -0.577350269189626;
		qp_coords(2, 0) = -0.577350269189626;
		qp_coords(2, 1) = 0.577350269189626;
		qp_coords(3, 0) = 0.577350269189626;
		qp_coords(3, 1) = 0.577350269189626;
		return qp_coords;
	}

	Eigen::Vector4d make_weights()
	{
		Eigen::Vector4d weight;
		weight[0] = 1.0;
		weight[1] = 1.0;
		weight[2] = 1.0;
		weight[3] = 1.0;
		return weight;
	}

	const Eigen::Matrix<double, 4, 2> CQuadrature::qp_coords = make_quadrature();
	const Eigen::Vector4d CQuadrature::weights = make_weights();

}

