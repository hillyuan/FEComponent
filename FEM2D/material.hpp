#ifndef ROLLFEM2D_MATERIAL_HPP
#define ROLLFEM2D_MATERIAL_HPP

#include <iostream>
#include <valarray>
#include <string>

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
		void StressUpdate(double strain[3], double stress[3]) const;
	};
}


#endif
