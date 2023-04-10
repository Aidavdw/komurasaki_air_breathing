#pragma once
#include <stdexcept>
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
		return { x + other.x, y+other.y };
	}

	inline Position operator * (const double scale) const
	{
		return {x*scale, y*scale};
	}

	inline Position operator - (const Position& other) const
	{
		return {x - other.x, y - other.y};
	}
};