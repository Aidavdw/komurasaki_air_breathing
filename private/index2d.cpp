#include "index2d.h"

#include <string>

#include "AuxFunctions.h"

std::string CellIndex::ToString() const
{
    return "(x=" + std::to_string(x) + ", y=" + std::to_string(y) + "), up=" + FaceToString(upDirection);
}

CellIndex TransformToOtherCoordinateSystem(const CellIndex& positionInOtherCoordinateSystem,
                                           const CellIndex& fromOrigin, const CellIndex& toOrigin)
{
    /*
    *				    TOP
    *		 + -- > --- > --- > --- +
    *	L	 |			|			|	R
    *	E	/\			\/			\/	I
    *	F	 | ->				 <- |	G
    *	T	/\          /\			\/	H
    *		 |          |			|	T
    *		 + -- < --- < --- < --- +
    *				  BOTTOM
    */	 
        // The order of operations is important here; First de-rotate the original coordinate system

    // Express the position as what it would be if the other coordinate system was oriented in the same direction as the target one.
    // In essence, this is doing an inverse rotation in local space.
    CellIndex posInLocalDeRotatedReferenceFrame;
     const int amountOfCounterClockwiseQuarterRotations = AmountOfNinetyDegreeRotationsBetweenOrientations(fromOrigin.upDirection, toOrigin.upDirection); // The amount of counterclockwise rotations that the original coordinate system has to do to align with the new one.

    // Now that the rotation required is set, apply it
    switch (amountOfCounterClockwiseQuarterRotations)
    {
    case 0:
        break;
    case 1:
        posInLocalDeRotatedReferenceFrame = {-positionInOtherCoordinateSystem.y, positionInOtherCoordinateSystem.x, toOrigin.upDirection};
        break;
    case 2:
        posInLocalDeRotatedReferenceFrame = {-positionInOtherCoordinateSystem.x, -positionInOtherCoordinateSystem.y, toOrigin.upDirection};
        break;
    case -1:
        posInLocalDeRotatedReferenceFrame = {positionInOtherCoordinateSystem.y, -positionInOtherCoordinateSystem.x, toOrigin.upDirection};
        break;
    default:
        throw std::logic_error("Impossible transformation in converting coordinate systems for 2 Pos objects.");
    }
    

    return toOrigin + posInLocalDeRotatedReferenceFrame;
}
