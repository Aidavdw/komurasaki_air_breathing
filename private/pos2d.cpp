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