#pragma once
#include <stdexcept>
#include "domain_enums.h"
// Really simple wrapper that represents a position in a domain
struct Position
{
	Position() :
		x(0),
		y(0),
		upDirection(EBoundaryLocation::TOP)
	{}

	Position(const double x, const double y, const EBoundaryLocation upDirection) :
		x(x),
		y(y),
		upDirection(upDirection)
	{}

	Position(const double x, const double y) :
	x(x),
	y(y),
	upDirection(EBoundaryLocation::TOP)
	{}

	double x;
	double y;

	EBoundaryLocation upDirection; // useful for when using rotating reference frames. Determines what is up.

	// Returns the distance between {0,0} and this position (when considering it as a vector).
	inline double Distance() const;

	inline double& operator [] (int axis)
	{
		if (axis == 0)
			return x;
		else if (axis == 1)
			return y;
		else
			throw std::invalid_argument("not a valid axis for a 2d position.");
	}

	
	inline Position operator + (const Position& other) const
	{
		#ifdef _DEBUG
		// If two positions are given in different coordinate frames, use PlusPositionInOtherCoordinateFrame(). Deliberately separated to make it easier to catch bugs.
		if (upDirection != other.upDirection)
			throw std::logic_error("Cannot add two positions in different coordinate frames using the + operator, you have to explicitly use PlusPositionInOtherCoordinateFrame(). This is to make it easier to catch bugs :)");
		#endif
		
		return { x + other.x, y+other.y, upDirection};
	}

	inline Position operator * (const double scale) const
	{
		return {x*scale, y*scale, upDirection};
	}

	inline Position operator - (const Position& other) const
	{
		return operator+(other * -1) ;
	}
};

Position TransformToOtherCoordinateSystem(const Position& positionInOtherCoordinateSystem, const Position& fromOrigin, const Position& toOrigin);