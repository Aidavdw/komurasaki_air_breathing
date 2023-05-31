#pragma once
#include <stdexcept>
#include "domain_enums.h"

// Really simple wrapper that represents a cell index in a grid or matrix.
struct CellIndex
{
	CellIndex() :
		x(0),
		y(0),
		relativeToBoundary(BOTTOM)
	{}

	
	CellIndex(const int x,const int y) :
		x(x),
		y(y),
		relativeToBoundary(BOTTOM)
	{}

	CellIndex(const int x,const int y, const EFace datumBoundary) :
		x(x),
		y(y),
		relativeToBoundary(datumBoundary)
	{}

	int x;
	int y;
	EFace relativeToBoundary; // useful for when using rotating reference frames. Determines what is up.

	// Check if actually used?
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
#ifdef _DEBUG
		// Ensure that they have the same up-direction
		if (relativeToBoundary != other.relativeToBoundary)
			throw std::logic_error("Cannot add two positions in different coordinate frames using the + operator, you have to explicitly use PlusPositionInOtherCoordinateFrame(). This is to make it easier to catch bugs :)");
#endif
		
		return { x + other.x, y+other.y };
	}

	std::string ToString() const;

};

CellIndex TransformToOtherCoordinateSystem(const CellIndex& positionInOtherCoordinateSystem, const CellIndex& fromOrigin, const CellIndex& toOrigin);
