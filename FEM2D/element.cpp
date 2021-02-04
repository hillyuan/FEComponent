#include "element.hpp"

#include <fstream>
#include <sstream>

namespace ROLLFEM2D
{
	void ShapeFunc(double lcoord[2], double sfunc[4])
	{
		sfunc[0] = 0.25 * (1.0 - lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[1] = 0.25 * (1.0 + lcoord[0]) * (1.0 - lcoord[1]);
		sfunc[2] = 0.25 * (1.0 + lcoord[0]) * (1.0 + lcoord[1]);
		sfunc[3] = 0.25 * (1.0 - lcoord[0]) * (1.0 + lcoord[1]);
	}

	void ShapeDeriv(double lcoord[2], double sfunc[4][2])
	{
		sfunc[0][0] = -0.25 * (1.0 - lcoord[1]);
		sfunc[1][0] = 0.25 * (1.0 - lcoord[1]);
		sfunc[2][0] = 0.25 * (1.0 + lcoord[1]);
		sfunc[3][0] = -0.25 * (1.0 + lcoord[1]);

		sfunc[0][1] = -0.25 * (1.0 - lcoord[0]);
		sfunc[1][1] = 0.25 * (1.0 - lcoord[0]);
		sfunc[2][1] = 0.25 * (1.0 + lcoord[0]);
		sfunc[3][1] = -0.25 * (1.0 + lcoord[0]);
	}

}

