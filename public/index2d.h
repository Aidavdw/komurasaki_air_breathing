#pragma once
#include <stdexcept>

// Really simple wrapper that represents a cell index in a grid or matrix.
struct CellIndex
{
	CellIndex(const int x,const int y) :
		x(x),
		y(y)
	{}

	int x;
	int y;

	inline int& operator [] (const int axis)
	{
		if (axis == 0)
			return x;
		else if (axis == 1)
			return y;
		else
			throw std::invalid_argument("Invalid axis");
	}

	inline CellIndex operator + (const CellIndex& other) const
	{
		return { x + other.x, y+other.y };
	}
};