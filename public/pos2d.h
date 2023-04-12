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

	double x;
	double y;

	EBoundaryLocation upDirection; // useful for when using rotating reference frames. Determines what is up.

	// Returns the distance between {0,0} and this position (when considering it as a vector).
	inline double Distance() const;

	// Adds two numbers in different coordinate frames. Returns the answer in the coordinate frame of the object that this is called from.
	Position PlusPositionInOtherCoordinateFrame(const Position& other) const;

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
		// If two positions are given in different coordinate frames, use PlusPositionInOtherCoordinateFrame(). Deliberately separated to make it easier to catch bugs.
		if (upDirection == other.upDirection)
			return { x + other.x, y+other.y, upDirection};
		else
			throw std::logic_error("Cannot add two positions in different coordinate frames using the + operator, you have to explicitly use PlusPositionInOtherCoordinateFrame(). This is to make it easier to catch bugs :)");
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