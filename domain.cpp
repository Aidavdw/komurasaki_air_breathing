#include "domain.h"
#include "sim_case.h"
#include <stdexcept>



void Domain::SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type)
{
	if (type == EBoundaryType::CONNECTED)
	{
		throw std::invalid_argument("Boundary type cannot be manually set to connected. Use ConnectBoundaries() instead.");
	}
	// Use the ordering of the enum to index the boundaries[] array
	boundaries[location] = Boundary(type);
}
