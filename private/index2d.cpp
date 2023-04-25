#include "index2d.h"

CellIndex TransformToOtherCoordinateSystem(const CellIndex& positionInOtherCoordinateSystem,
    const CellIndex& fromOrigin, const CellIndex& toOrigin)
{
        // The order of operations is important here; First de-rotate the original coordinate system
    
     int amountOfCounterClockwiseQuarterRotations; // The amount of counterclockwise rotations that the original coordinate system has to do to align with the new one.

    // Not the most elegant way to do this, but doesn't make any assumptions on the layout of the struct in its cpp definition.
    switch (fromOrigin.upDirection)
    {
    case TOP:
        switch (toOrigin.upDirection)
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
        switch (toOrigin.upDirection)
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
        switch (toOrigin.upDirection)
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
        switch (toOrigin.upDirection)
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

    // Express the position as what it would be if the other coordinate system was oriented in the same direction as the target one.
    // In essence, this is doing an inverse rotation in local space.
    CellIndex posInLocalDeRotatedReferenceFrame;

    // Now that the rotation required is set, apply it
    switch (amountOfCounterClockwiseQuarterRotations)
    {
    case 0:
        break;
    case 1:
        posInLocalDeRotatedReferenceFrame = {positionInOtherCoordinateSystem.y, positionInOtherCoordinateSystem.x};
        break;
    case 2:
        posInLocalDeRotatedReferenceFrame = {-positionInOtherCoordinateSystem.x, -positionInOtherCoordinateSystem.y};
        break;
    case -1:
        posInLocalDeRotatedReferenceFrame = {-positionInOtherCoordinateSystem.y, -positionInOtherCoordinateSystem.x};
        break;
    default:
        throw std::logic_error("Impossible transformation in converting coordinate systems for 2 Pos objects.");
    }

    return fromOrigin + posInLocalDeRotatedReferenceFrame;
}
