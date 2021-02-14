#ifndef ROLLFEM2D_CONTROL_HPP
#define ROLLFEM2D_CONTROL_HPP

#include <map>

#include "mesh.hpp"
#include "material.hpp"
#include "boundary.hpp"

namespace ROLLFEM2D
{

	class CControl
	{
	private:
		CMesh mesh;
		Eigen::SparseMatrix<double> StiffMatrix;
		Eigen::VectorXd				loads;

		std::vector< Constraint > constraints;
		std::vector< CLoad > cloads;
		std::vector< DLoad > dloads;
		
	public:
		CControl::CControl(char* file);

		void calGlobalStiffMatrix()
		{
			mesh.calGlobalStiffMatrix(StiffMatrix);
		}
		void ApplyConstraints();
		void ApplyDistributedLoads();
		void solve();

	};
}


#endif
