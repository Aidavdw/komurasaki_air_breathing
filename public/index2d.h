#pragma once
#include <stdexcept>
#include "domain_enums.h"

// Really simple wrapper that represents a cell index in a grid or matrix.
struct CellIndex
{
	CellIndex() :
		x(NAN),
		y(NAN),
		upDirection(TOP)
	{}

	
	CellIndex(const int x,const int y) :
		x(x),
		y(y),
		upDirection(TOP)
	{}

	CellIndex(const int x,const int y, const EBoundaryLocation upDirection) :
		x(x),
		y(y),
		upDirection(upDirection)
	{}

	int x;
	int y;
	EBoundaryLocation upDirection; // useful for when using rotating reference frames. Determines what is up.

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

CellIndex TransformToOtherCoordinateSystem(const CellIndex& positionInOtherCoordinateSystem, const CellIndex& fromOrigin, const CellIndex& toOrigin);
