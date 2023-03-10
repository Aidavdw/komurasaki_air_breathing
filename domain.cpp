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

void Domain::InitialiseDomain()
{
	NXtot = std::max(2 + 2 * NGHOST, (XLENGTH[i] / dx / graded_ratio_x[i] + 2 * NGHOST));
	NYtot = std::max(2 + 2 * NGHOST, (YLENGTH[i] / dy / graded_ratio_y[i] + 2 * NGHOST));
}
