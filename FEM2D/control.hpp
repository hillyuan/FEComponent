#ifndef ROLLFEM2D_CONTROL_HPP
#define ROLLFEM2D_CONTROL_HPP

#include "mesh.hpp"
#include "material.hpp"

namespace ROLLFEM2D
{
	class CControl
	{
	private:
		CMesh mesh;
		Eigen::SparseMatrix<double> StiffMatrix;
		
	public:
		CControl::CControl(char* file);

		void calGlobalStiffMatrix()
		{
			mesh.calGlobalStiffMatrix(StiffMatrix);
		}

	};
}


#endif
