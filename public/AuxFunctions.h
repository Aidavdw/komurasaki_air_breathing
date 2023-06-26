#pragma once
/* This file contains several auxiliary functions that are used throughout the code.
 */
#include <limits>
#include <utility>
#include <vector>

#include "pos2d.h"

# define M_PI 3.14159265358979323846


bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() ); // doing (double == double) is dangerous because of machine epsilon. This gets a 'good enough' estimate when used in combination with subtraction.

std::pair<Position, Position> ExtrudeAlongNormal(const Position startPos, const Position endPos, const double depth); // Takes a line between two points, and extends it orthogonally up to a depth, creating a square. The output pair are the two newly created points.

size_t FindIndexLeftOfValueByBisection(const std::vector<double>& field, const double valueToFind); // Given a vector of doubles, returns the index of the element  that is the closest to- while also being lower than the value to find. [10,11,12,13,14], to find = 12.5 -> returns 2

bool IsCloserToLeftThanToRight(const double valueToFind, const double lValue, const double rValue, const bool bOnlyAllowBetweenTheseTwoValues=false); // Given two reference values, checks which one it is closest to. input: 20, l=15, r=40 -> return true;

int AmountOfNinetyDegreeRotationsBetweenOrientations(const EFace from, const EFace to); // Returns how many times you would have to rotate the 'from' direction to get to 'to'.

std::pair<double, double> GetBoundsWithDelta(const std::vector<double>& sortedVector); // Returns the inside-extrapolated left- and right bound of the vector.