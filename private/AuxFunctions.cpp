﻿#include "AuxFunctions.h"

#include <algorithm>
#include <cassert>
#include <string>

#include "pos2d.h"

#define MATCHING_TOLERANCE 1.0E-4


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
    auto maxOfAll = std::max_element(field.begin(), field.end());
    auto minOfAll = std::min_element(field.begin(), field.end());
    if (!IsCloseToZero(*maxOfAll - field.back()) || !IsCloseToZero(*minOfAll - field.front()))
        throw std::logic_error("FindIndexLeftOfValueByBisection() can only be used on an array that is increasing left to right.");
    if (!std::is_sorted(field.begin(), field.end()))
        throw std::logic_error("FindIndexLeftOfValueByBisection() can only be used on an array that is sorted.");
#endif

    double min = field.front();
    double max = field.back();
    size_t l = 0;   // The current left index.
    size_t r = field.size() - 1;   // The current right index.   

    // Ensure that the value to find is actually between the lValue and the rValue. This assumes the field is ever increasing (in either direction), and hence the values at the borders are the highest.
    // If it's just slightly above or below the max/min, also allow it.
    auto withDelta = GetBoundsWithDelta(field);
    if (valueToFind < withDelta.first || valueToFind > withDelta.second)
        throw std::range_error("Value to find should be between the max and min value of the vector ( " + std::to_string(max) +"-> " + std::to_string(min) + ")");

#ifdef _DEBUG
    int nIters = 0; // Debug value to make sure the bisection does not loop indefinitely.
#endif
    
    while (r-l != 1)
    {
#ifdef _DEBUG
        if (nIters > 1000)
            throw std::range_error("Bisection took more than 1000 tries to converge. Are you doing this right?");
        nIters++;
        // todo: implement test to see if the array is actually sorted. > or < does not flip
#endif

        size_t c = l + static_cast<size_t>((r - l)/2); // The index of the cell in between l and r
        double cValue = field.at(c);
        if (valueToFind < cValue)
            r = c;
        else
            l = c;
    }
        return l;
}

bool IsCloserToLeftThanToRight(const double valueToFind, const double lValue, const double rValue, const bool bOnlyAllowBetweenTheseTwoValues)
{
    if (bOnlyAllowBetweenTheseTwoValues)
    {
        // Ensure that the value to find is actually between the lValue and the rValue. This assumes the field is ever increasing (in either direction), and hence the values at the borders are the highest. Applies a delta to deal with floating-point.
        if (valueToFind < (lValue - MATCHING_TOLERANCE) || valueToFind > (rValue + MATCHING_TOLERANCE))
            throw std::range_error("Value to find is outside of the extremes of the input vector.");
    }
    
    // See which of two points it is closer to
    double avg = (lValue + rValue)/2;
    if (valueToFind < avg)
        return true;
    else
        return false;
}

int AmountOfNinetyDegreeRotationsBetweenOrientations(const EFace from, const EFace to)
{
    int amountOfCounterClockwiseQuarterRotations;
        // Not the most elegant way to do this, but doesn't make any assumptions on the layout of the struct in its cpp definition.
    switch (from)
    {
    case TOP: // original up
        switch (to)
        {
            case TOP: // new one: right
                amountOfCounterClockwiseQuarterRotations = 0;
                break;
            case LEFT:
                amountOfCounterClockwiseQuarterRotations = 1;
                break;
            case RIGHT:
                amountOfCounterClockwiseQuarterRotations = -1;
                break;
            case BOTTOM:
                amountOfCounterClockwiseQuarterRotations = 2;
                break;
        }
        break;
    case LEFT:
        switch (to)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        }
        break;
    case RIGHT:
        switch (to)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        }
        break;
    case BOTTOM:
        switch (to)
        {
        case TOP:
            amountOfCounterClockwiseQuarterRotations = 2;
                break;
        case LEFT:
            amountOfCounterClockwiseQuarterRotations = -1;
                break;
        case RIGHT:
            amountOfCounterClockwiseQuarterRotations = 1;
                break;
        case BOTTOM:
            amountOfCounterClockwiseQuarterRotations = 0;
                break;
        }
        break;
    }

#ifdef _DEBUG
    if (amountOfCounterClockwiseQuarterRotations < -1 || amountOfCounterClockwiseQuarterRotations > 2)
        throw std::logic_error("Value uninitialised in AmountOfNinetyDegreeRotationsBetweenOrientations()");
#endif

    return amountOfCounterClockwiseQuarterRotations;
}

std::pair<double, double> GetBoundsWithDelta(const std::vector<double>& field)
{
    double minWithDelta = field.front() - std::abs(field.front() + field.at(1))/2 - MATCHING_TOLERANCE;
    double maxWithDelta = field.back() + std::abs(field.back() + field.at(field.size()-2) )/2 + MATCHING_TOLERANCE;

    return {minWithDelta, maxWithDelta};
}

EAxisDirection GetUpDirectionOfBoundary(const EFace boundary)
{
    if (boundary == TOP || boundary == BOTTOM)
        return EAxisDirection::POSITIVE;
    if (boundary == LEFT || boundary == RIGHT)
        return EAxisDirection::NEGATIVE;
    else
        throw std::invalid_argument("Cannot give up direction of this type of boundary");
}

void SwapValues(double& a, double& b)
{
    const double buf = a;
    a = b;
    b = buf;
}
