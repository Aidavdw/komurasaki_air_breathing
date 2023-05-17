#include "AuxFunctions.h"
#include <cassert>


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

size_t FindIndexLeftOfValueByBisection(const std::vector<double>& field, const double valueToFind)
{
#ifdef _DEBUG
    assert(!field.empty());
#endif

    size_t l = 0;   // The current left index.
    size_t r = field.size() - 1;   // The current right index.
    
    double lValue = field.front();
    double rValue = field.back();

    bool bLeftLowerThanRight;
    if (lValue < rValue)
    {
        bLeftLowerThanRight = true;
        // Ensure that the value to find is actually between the lValue and the rValue. This assumes the field is ever increasing (in either direction), and hence the values at the borders are the highest.
        if (valueToFind < lValue || valueToFind > rValue)
            throw std::range_error("Value to find is outside of the extremes of the input vector.");
    }
    else // lValue < rValue
    {
        bLeftLowerThanRight = false;
        if (valueToFind < rValue || valueToFind > lValue)
            throw std::range_error("Value to find is outside of the extremes of the input vector.");
    }

#ifdef _DEBUG
    int nIters = 0; // Debug value to make sure the bisection does not loop indefinitely.
#endif
    
    while (r-l != 1)
    {
#ifdef _DEBUG
        if (nIters > 1000)
            throw std::range_error("Bisection took more than 1000 tries to converge. Are you doing this right?");
        nIters++;
#endif

        size_t c = l + static_cast<size_t>((r - l)/2); // The index of the cell in between l and r
        double cValue = field.at(c);

        if (bLeftLowerThanRight) // lValue < rValue
        {
            if (valueToFind < cValue)
                r = c;
            else
                l = c;
        }
        else // lValue > rValue
        {
            if (valueToFind > cValue)
                r = c;
            else
                l = c;
        }
    }
        return l;

}

bool IsCloserToLeftThanToRight(const double valueToFind, const double lValue, const double rValue)
{
    bool bLeftLowerThanRight;
    if (lValue < rValue)
    {
        bLeftLowerThanRight = true;
        // Ensure that the value to find is actually between the lValue and the rValue. This assumes the field is ever increasing (in either direction), and hence the values at the borders are the highest.
        if (valueToFind < lValue || valueToFind > rValue)
            throw std::range_error("Value to find is outside of the extremes of the input vector.");
    }
    else // lValue < rValue
        {
        bLeftLowerThanRight = false;
        if (valueToFind < rValue || valueToFind > lValue)
            throw std::range_error("Value to find is outside of the extremes of the input vector.");
        }
    
    // See which of two points it is closer to
    double avg = (lValue + rValue)/2;
    if (bLeftLowerThanRight) // lValue < rValue
        {
        if (valueToFind < avg)
            return true;
        else
            return false;
        }
    else // lValue > rValue
        {
        if (valueToFind > avg)
            return true;
        else
            return false;
        }
}
