#pragma once
#include <limits>
#include <utility>

#include "pos2d.h"


bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() );

std::pair<Position, Position> ExtrudeAlongNormal(const Position startPos, const Position endPos, const double depth);

