// A small struct containing all the quantities in one cell. Differently organised from that in the domain, but this is a good way to transfer some data more easily.
struct CellValues
{
	CellValues() :
		rho(0),
		u(0),
		v(0),
		p(0),
		E(0),
		T(0),
		H(0)
	{};

	double rho; // Density
	double u; // velocity x-component
	double v; // velocity y-component
	double p; // pressure
	double E; // Internal energy?
	double T; // Temperature
	double H; // enthalpy ?
};