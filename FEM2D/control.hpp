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

		std::vector< Constraint > constraints;
		
	public:
		CControl::CControl(char* file);

		void calGlobalStiffMatrix()
		{
			mesh.calGlobalStiffMatrix(StiffMatrix);
		}
		void ApplyConstraints();

	};
}


#endif
