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
		std::string name;
		double lame_lambda;
		double lame_mu;
		double Youngs;
		double Poission;
		
	public:
		CMaterial(std::string&, double&, double&);

		Eigen::Matrix<double, 3, 3> ElasticMatrix;
		void StressUpdate(double strain[3], double stress[3]) const;

		void print(std::ostream& os) const;
	};
}


#endif
