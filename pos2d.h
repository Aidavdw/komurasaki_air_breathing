
// Really simple wrapper that represents a position in a domain
struct PositionInDomain
{
	PositionInDomain(const double x, const double y) :
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