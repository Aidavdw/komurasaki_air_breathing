#include "AuxFunctions.h"


bool IsCloseToZero(const double x, const double tolerance)
{
    return std::abs(x) < tolerance;
}

std::pair<Position, Position> ExtrudeAlongNormal(const Position startPos, const Position endPos, const double depth)
{
    /*
     *       endPos + ============= + extrudedEndPos
     *              |               ||
     *              |               ||
     *    startPos  + ============= + extrudedStartPos
     */

    const Position edge = (endPos - startPos);
    Position normal = Position(edge.y, edge.x) * (1/edge.Distance());
    Position extrudedStartPos = startPos + normal * depth;
    Position extrudedEndPosEndPos = extrudedStartPos + edge;

    return std::make_pair(extrudedStartPos, extrudedEndPosEndPos);
}
