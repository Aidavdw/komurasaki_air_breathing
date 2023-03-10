#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <algorithm>
#include "parameters.h"


void Domain::SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type)
{
	if (type == EBoundaryType::CONNECTED)
	{
		throw std::invalid_argument("Boundary type cannot be manually set to connected. Use ConnectBoundaries() instead.");
	}
	// Use the ordering of the enum to index the boundaries[] array
	boundaries[location] = Boundary(type);
}

int Domain::GetAmountOfCellsInAxis(const unsigned int axis) const
{
	ValidateAxisInput(axis);
	std::max(2 + 2 * NGHOST, gridResolution[axis] + 2 * NGHOST);
}

int Domain::GetCellResolutionInAxis(const int axis) const
{
	ValidateAxisInput(axis);
	return size[axis] / gridResolution[axis];
}

int Domain::GetTotalAmountOfCells() const
{
	return gridResolution[0] * gridResolution[y]
}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y!");
	}
}
