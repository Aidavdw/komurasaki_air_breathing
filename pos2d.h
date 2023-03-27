
// Really simple wrapper that represents a position in a domain
struct Position
{
	Position() :
		x(0),
		y(0)
	{};

	Position(const double x, const double y) :
		x(x),
		y(y)
	{};

	double x;
	double y;

	inline double& operator [] (int axis)
	{
		if (axis == 0)
			return x;
		else if (axis == 1)
			return y;
	}
};