#include "AuxFunctions.h"


bool IsCloseToZero(const double x, const double tolerance)
{
    return std::abs(x) < tolerance;
}