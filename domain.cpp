#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <algorithm>
#include "parameters.h"


Domain::Domain(std::string& name, const double position_arg[2], const double size_arg[2], const int amountOfCells_arg[2]) :
	name(name)
{
	// Not a very pretty way to do this, but initialiser lists appear to break when using c style arrays
	position[0] = position_arg[0];
	position[1] = position_arg[1];
	size[0] = size_arg[0];
	size[1] = size_arg[1];
	amountOfCells[0] = amountOfCells_arg[0];
	amountOfCells[1] = amountOfCells_arg[1];

	// Initialising the field quantities with the given dimension
	rho = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	u = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	v = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	p = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	E = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	T = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	H = FieldQuantity(amountOfCells[0], amountOfCells[1]);

}

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
	std::max(2 + 2 * NGHOST, amountOfCells[axis] + 2 * NGHOST);
}

int Domain::GetCellResolutionInAxis(const int axis) const
{
	ValidateAxisInput(axis);
	return size[axis] / amountOfCells[axis];
}

int Domain::GetTotalAmountOfCells() const
{
	return amountOfCells[0] * amountOfCells[1];
}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y!");
	}
}
