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
		int type;              // 0: volume force; 1: surface pressuure 2: surface load in direction[2]
		double val;            // surface pressure; maybe a table
		double direction[2];

		DLoad(std::string& name, std::vector<double>& v)
			: SetName(name)
		{
			if (v.size() == 1) {
				val = v[0]; type = 1;
			}
			else {
				direction[0] = v[0]; direction[1] = v[1];
				val = 1.0;
				type = 2;
			}
		};
	};
}


#endif
