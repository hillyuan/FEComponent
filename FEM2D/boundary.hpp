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
		std::string NSetName;   // node set name
		int id_node;    // nodal index
		Type type;
		double val;

		Constraint(std::string& name, int& t, double& v)
			: NSetName(name), val(v), id_node(-1)
		{
			type = static_cast<Type>(t);
		};
		Constraint(int& id, int& t, double& v)
			: id_node(id), val(v)
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
		std::string SetName;      // edge/element sets name
		int type;                 // 0: volume force; 1: surface pressuure 2: surface load in direction[2]
		double val;               // surface pressure
		std::vector<double> pos;  
		std::vector<double> vals; // surface pressure upon pos, type must be 1 in this case
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

		DLoad(std::string& name) : SetName(name)
		{};
	};

	// initial strain
	struct InitStrain
	{
		std::string SetName;      // element set name
		Eigen::Vector3d strain;   // strain value

		InitStrain(const std::string& name, double s0, double s1, double s2)
			: SetName(name)
		{
			strain = Eigen::Vector3d(s0, s1, s2);
		}
	};
}


#endif
