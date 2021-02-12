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
		std::vector<CMaterial> materials;
		
	public:
		CControl::CControl(std::string& file);

	};
}


#endif
