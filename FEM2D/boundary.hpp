#ifndef ROLLFEM2D_BOUNDARY_HPP
#define ROLLFEM2D_BOUNDARY_HPP


namespace ROLLFEM2D
{
	// Displacement constraint
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

	// Concerntrated loads
	struct CLoad
	{
		std::string NSetName;   // node set name
		double fx,fy;           // force value
	};

	// Distributed loads
	struct DLoad
	{
		std::string SetName;   // edge/element sets name
		int type;              // 0: volume force; 1: surface pressuure
		double val;            // surface pressure; maybe a table

		DLoad(std::string& name, double& v)
			: SetName(name), val(v), type(0)
		{};
	};
}


#endif
