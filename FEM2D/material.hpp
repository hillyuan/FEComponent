#ifndef ROLLFEM2D_MATERIAL_HPP
#define ROLLFEM2D_MATERIAL_HPP

#include <iostream>
#include <valarray>
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

		double ElasticMatrix[3][3];
		void StressUpdate(std::valarray<double>&, std::valarray<double>&) const;
	};
}


#endif
