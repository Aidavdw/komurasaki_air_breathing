#include "index2d.h"

#include <string>

#include "AuxFunctions.h"

std::string CellIndex::ToString() const
{
    return "(x=" + std::to_string(x) + ", y=" + std::to_string(y) + "), up=" + FaceToString(relativeToBoundary);
}

CellIndex TransformToOtherCoordinateSystem(const CellIndex& positionInOtherCoordinateSystem,
                                           const CellIndex& fromOrigin, const CellIndex& toOrigin)
{
    /*
    *				    TOP
    *					/\
    *					|
    *		 + -- > --- > --- > --- +
    *	L	 |			 			|		R
    *	E	/\			 			/\		I
    *	F	 | ->					|  ->	G
    *	T	/\          /\			/\		H
    *		 |          |			|		T
    *		 + -- > --- > --- > --- +
    *				  BOTTOM
    */
        // The order of operations is important here; First de-rotate the original coordinate system

#ifdef _DEBUG
    if (fromOrigin.relativeToBoundary != positionInOtherCoordinateSystem.relativeToBoundary)
        throw std::logic_error("positionInOtherCoordinate system must have the same orientation as FromOrigin!");
#endif
    
    const bool bSameAxis = (fromOrigin.relativeToBoundary == toOrigin.relativeToBoundary || fromOrigin.relativeToBoundary == Opposite(fromOrigin.relativeToBoundary));
    if (bSameAxis) 
    {
        // If their up-axis are aligned, only the origin/datum is different.
        return {positionInOtherCoordinateSystem.x - fromOrigin.x, positionInOtherCoordinateSystem.y - fromOrigin.y, toOrigin.relativeToBoundary};
    }
    else // It's pointing in the other direction. Since both LEFT and RIGHT point to the positive axis (->) and TOP and BOTTOM point to the positive axis ( /\ ), any x is just y, and any y = x. No weird minuses, as the local frame is left-handed.
    {
        return {positionInOtherCoordinateSystem.y - fromOrigin.x, positionInOtherCoordinateSystem.x - fromOrigin.y, toOrigin.relativeToBoundary};
    }
}
