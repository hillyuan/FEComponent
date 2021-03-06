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
		int type;             // 0: plane stress; 1: plane strain
		double lame_lambda;
		double lame_mu;
		double Youngs;
		double Poission;
		double Density;
		
	public:
		CMaterial(int pbtype, std::string&, double&, double&, double&);

		Eigen::Matrix<double, 3, 3> ElasticMatrix;
		Eigen::Vector3d StressUpdate(Eigen::Vector3d&) const;

		void print(std::ostream& os) const;
		double getDensity() const
		{
			return Density;
		}
	};
}


#endif
