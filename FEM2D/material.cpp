#include "material.hpp"

#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	CMaterial::CMaterial(int pbtype, std::string& name, double& y, double& p)
		: name(name), Youngs(y), Poission(p), type(pbtype), density(0.0)
	{
		lame_lambda = y*p/(1+p)/(1-2.0*p);
		lame_mu = 0.5*y/(1+p);

		if (type == 0)
		{
			double cc = Youngs / (1.0 - Poission * Poission);
			ElasticMatrix(0, 0) = cc;
			ElasticMatrix(0, 1) = cc * Poission;
			ElasticMatrix(0, 2) = 0.0;
			ElasticMatrix(1, 0) = ElasticMatrix(0, 1);
			ElasticMatrix(1, 1) = ElasticMatrix(0, 0);
			ElasticMatrix(1, 2) = 0.0;
			ElasticMatrix(2, 0) = 0.0;
			ElasticMatrix(2, 1) = 0.0;
			ElasticMatrix(2, 2) = 0.5 * Youngs / (1.0 + Poission);
		}
		else {
			ElasticMatrix(0, 0) = lame_lambda + 2.0 * lame_mu;
			ElasticMatrix(0, 1) = lame_lambda;
			ElasticMatrix(0, 2) = 0.0;
			ElasticMatrix(1, 0) = lame_lambda;
			ElasticMatrix(1, 1) = ElasticMatrix(0, 0);
			ElasticMatrix(1, 2) = 0.0;
			ElasticMatrix(2, 0) = 0.0;
			ElasticMatrix(2, 1) = 0.0;
			ElasticMatrix(2, 2) = lame_mu;
		}
	}

	Eigen::Vector3d CMaterial::StressUpdate(Eigen::Vector3d& strain) const
	{
		Eigen::Vector3d stress = ElasticMatrix * strain;

	//	stress[0] = 2.0 * lame_mu * strain[0] + lame_lambda * (strain[0] + strain[1]);
	//	stress[1] = 2.0 * lame_mu * strain[1] + lame_lambda * (strain[0] + strain[1]);
	//	stress[2] = lame_mu * strain[2];

		return stress;
	}

	void CMaterial::print(std::ostream& os) const
	{
		std::cout << "Material: " << name << std::endl;
		std::cout << "  Youngs modulus= " << Youngs << " ; Poisson's ratio=" << Poission << std::endl;
	}

}

