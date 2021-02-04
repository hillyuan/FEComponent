#include "material.hpp"

#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	CMaterial::CMaterial(double& y, double& p)
		: Youngs(y), Poission(p)
	{
		lame_lambda = y*p/(1+p)/(1-2.0*p);
		lame_mu = 0.5*y/(1+p);

		ElasticMatrix[0][0] = lame_lambda + 2.0 * lame_mu;
		ElasticMatrix[0][1] = lame_lambda;
		ElasticMatrix[0][2] = 0.0;
		ElasticMatrix[1][0] = lame_lambda;
		ElasticMatrix[1][1] = ElasticMatrix[0][0];
		ElasticMatrix[1][2] = 0.0;
		ElasticMatrix[2][0] = 0.0;
		ElasticMatrix[2][1] = 0.0;
		ElasticMatrix[2][2] = lame_mu;
	}

	void CMaterial::StressUpdate(std::valarray<double>& strain, std::valarray<double>& stress) const
	{
		stress[0] = 2.0 * lame_mu * strain[0] + lame_lambda * (strain[0] + strain[1]);
		stress[1] = 2.0 * lame_mu * strain[1] + lame_lambda * (strain[0] + strain[1]);
		stress[2] = lame_mu * strain[2];
	}

}
