#include "pos2d.h"
#include <cmath>

double Position::Distance() const
{
    // shorthand so that if either is zero, the square root will not be done to save performance.
    if (std::abs(x) < std::numeric_limits<double>::epsilon())
        return y;
    else if (std::abs(y) < std::numeric_limits<double>::epsilon())
        return x;

    return sqrt(x*x + y*y);
}

Position Position::PlusPositionInOtherCoordinateFrame(const Position& other) const
{
    int amountOfCounterClockwiseQuarterRotations; // The amount of counterclockwise rotations that the other object needs to do to match the master position.

    // Not the most elegant way to do this, but doesn't make any assumptions on the layout of the struct in its cpp definition.
    switch (upDirection)
    {
    case TOP:
        switch (other.upDirection)
        {
            case TOP:
                amountOfCounterClockwiseQuarterRotations = 0;
                break;
            case LEFT:
                amountOfCounterClockwiseQuarterRotations = -1;
                break;
            case RIGHT:
                amountOfCounterClockwiseQuarterRotations = 1;
                break;
            case BOTTOM:
                amountOfCounterClockwiseQuarterRotations = 2;
                break;
        }
        break;
    case LEFT:
        switch (other.upDirection)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        }
        break;
    case RIGHT:
        switch (other.upDirection)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        }
        break;
    case BOTTOM:
        switch (other.upDirection)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        }
        break;
    }

    // Now that the rotation required is set, apply it
    switch (amountOfCounterClockwiseQuarterRotations)
    {
    case 0:
        return operator+(other);
    case 1:
        return {x + other.y, y + other.x, upDirection};
    case 2:
        return {x - other.x, y - other.y, upDirection};
    case -1:
        return {x - other.y, y - other.x, upDirection};
    default:
        throw std::logic_error("Impossible transformation in converting coordinate systems for 2 Pos objects.");
    }
}