
// Really simple wrapper that represents a cell index in a grid or matrix.
struct CellIndex
{
	CellIndex(const int x, int double y) :
		x(x),
		y(y)
	{};

	int x;
	int y;

	inline int& operator [] (int axis)
	{
		if (axis == 0)
			return x;
		else if (axis == 1)
			return y;
	}
};