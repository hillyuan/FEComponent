#ifndef ROLLFEM2D_CONTROL_HPP
#define ROLLFEM2D_CONTROL_HPP

#include <fstream>

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
		std::vector< InitStrain > initstrain;

		Eigen::VectorXd displacements;
		Eigen::VectorXd forces;

		std::string outfile;
		
	public:
		CControl::CControl(char* file);

		void calGlobalStiffMatrix()
		{
			mesh.calGlobalStiffMatrix(StiffMatrix);
		}
		void ApplyConstraints();
		void ApplyDistributedLoads();
		void ApplyInitialStrain();
		void Solve();
		void Update()
		{
			mesh.updateElements(displacements);
		}

		int VTKOutput() const;

	};
}


#endif
