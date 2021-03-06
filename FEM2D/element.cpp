#include "element.hpp"


namespace ROLLFEM2D
{
	Eigen::Vector<double,4> ShapeFunc(Eigen::Vector2d& lcoord)
	{
		Eigen::Vector<double,4> sfunc;
		sfunc[0] = 0.25 * (1.0 - lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[1] = 0.25 * (1.0 + lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[2] = 0.25 * (1.0 + lcoord[0]) * (1.0 + lcoord[1]);
		sfunc[3] = 0.25 * (1.0 - lcoord[0]) * (1.0 + lcoord[1]);
		return sfunc;
	}

	Eigen::Matrix<double, 4, 2> ShapeDeriv(Eigen::Vector2d& lcoord)
	{
		Eigen::Matrix<double, 4, 2> spderiv;
		spderiv(0,0) = -0.25 * (1.0 - lcoord[1]);
		spderiv(1,0) = 0.25 * (1.0 - lcoord[1]);
		spderiv(2,0) = 0.25 * (1.0 + lcoord[1]);
		spderiv(3,0) = -0.25 * (1.0 + lcoord[1]);

		spderiv(0,1) = -0.25 * (1.0 - lcoord[0]);
		spderiv(1,1) = -0.25 * (1.0 + lcoord[0]);
		spderiv(2,1) = 0.25 * (1.0 + lcoord[0]);
		spderiv(3,1) = 0.25 * (1.0 - lcoord[0]);
		return spderiv;
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

	std::array< Eigen::Vector<double, 4>, 4 > make_ShapeFuncs()
	{
		std::array< Eigen::Vector<double, 4>, 4 > spfunc;
		Eigen::Vector2d lcoord;
		for (unsigned int i = 0; i < 4; i++)
		{
			lcoord[0] = CQuadrature::qp_coords(i, 0);
			lcoord[1] = CQuadrature::qp_coords(i, 1);
			spfunc[i] = ShapeFunc(lcoord);
		}
		return spfunc;
	}

	std::array<Eigen::Matrix<double, 4, 2>, 4> make_ShapeDerivs()
	{
		std::array<Eigen::Matrix<double, 4, 2>, 4> spderivs;
		Eigen::Vector2d lcoord;
		for (unsigned int i = 0; i < 4; i++)
		{
			lcoord[0] = CQuadrature::qp_coords(i, 0);
			lcoord[1] = CQuadrature::qp_coords(i, 1);
			spderivs[i] = ShapeDeriv(lcoord);
		}
		return spderivs;
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

	Eigen::Vector<double, 8> CElement::calNodalForce() const
	{
		Eigen::Vector<double, 8> force;
		force.setZero();

		for (std::size_t npg = 0; npg < 4; npg++)
		{
			force += wg[npg] * stress[npg].transpose() * B[npg];
		}
		return force;
	}

	const Eigen::Matrix<double, 4, 2> CQuadrature::qp_coords = make_quadrature();
	const Eigen::Vector4d CQuadrature::weights = make_weights();
	const std::array<Eigen::Matrix<double, 4, 2>, 4> CQuadrature::ShapeDerivs = make_ShapeDerivs();
	const std::array< Eigen::Vector<double,4>, 4 > CQuadrature::ShapeFuncs = make_ShapeFuncs();
}

