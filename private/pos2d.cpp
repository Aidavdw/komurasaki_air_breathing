#include "pos2d.h"
#include <cmath>
#include <limits>

double Position::Distance() const
{
    // shorthand so that if either is zero, the square root will not be done to save performance.
    if (std::abs(x) < std::numeric_limits<double>::epsilon())
        return y;
    else if (std::abs(y) < std::numeric_limits<double>::epsilon())
        return x;

    return std::sqrt(x*x + y*y);
}


Position TransformToOtherCoordinateSystem(const Position& positionInOtherCoordinateSystem, const Position& fromOrigin, const Position& toOrigin)
{
    if (positionInOtherCoordinateSystem.upDirection != fromOrigin.upDirection)
        throw std::logic_error("position in other coordinate system and form origin must have the same up direction.");

    if (fromOrigin.upDirection == toOrigin.upDirection)
        return positionInOtherCoordinateSystem + (toOrigin - fromOrigin);
    else
    {
        // It is essentially flipped
        Position flippedPos = {positionInOtherCoordinateSystem.y, positionInOtherCoordinateSystem.x, toOrigin.upDirection};
        return {positionInOtherCoordinateSystem.y + (toOrigin.x - fromOrigin.x), positionInOtherCoordinateSystem.x + (toOrigin.y - fromOrigin.y), toOrigin.upDirection};
    }
}