#ifndef _ROLLFEM2D_ELEMENT_HPP
#define _ROLLFEM2D_ELEMENT_HPP

#include <iostream>
#include <valarray>
#include <string>
//#include <Eigen\Eigen>

namespace ROLLFEM2D
{
	void ShapeFunc(double lcoord[2], double sfunc[4]);

	void ShapeDeriv(double lcoord[2], double sfunc[4][2]);

	class CElement
	{
	public:
		int matl_id;
		std::size_t n0, n1, n2, n3;

	private:
		double shapefunc[4];
		double shapederiv[4][2];
	};
}


#endif
