#ifndef ROLLFEM2D_MATERIAL_HPP
#define ROLLFEM2D_MATERIAL_HPP

#include <iostream>
#include <string>
#include <Eigen\Eigen>

namespace ROLLFEM2D
{
	class CMaterial
	{
	private:
		double lame_lambda;
		double lame_mu;
		double Youngs;
		double Poission;
		
	public:
		CMaterial(double&, double&);

		Eigen::Matrix<double, 3, 3> ElasticMatrix;
		void StressUpdate(double strain[3], double stress[3]) const;
	};
}


#endif
