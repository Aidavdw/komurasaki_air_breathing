#pragma once
#include <limits>
#include <utility>
#include <vector>

#include "pos2d.h"

# define M_PI 3.14159265358979323846


bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() );

std::pair<Position, Position> ExtrudeAlongNormal(const Position startPos, const Position endPos, const double depth);

std::pair<size_t, size_t> FindIndexOfValueByBisection(const std::vector<double>& field, const double valueToFind);
