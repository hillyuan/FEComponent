#ifndef ROLLFEM2D_BOUNDARY_HPP
#define ROLLFEM2D_BOUNDARY_HPP


namespace ROLLFEM2D
{
	struct Constraint
	{
		enum Type
		{
			UX = 1 << 0,
			UY = 1 << 1,
			UXY = UX | UY
		};
		std::string NSetName;
		Type type;
		double val;

		Constraint(std::string& name, int& t, double& v)
			: NSetName(name), val(v)
		{
			type = static_cast<Type>(t);
		};
	};
}


#endif
