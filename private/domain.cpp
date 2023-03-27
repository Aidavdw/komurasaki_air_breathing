#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <algorithm>
#include "parameters.h"


Domain::Domain(const std::string& name, const Position position, const double sizeArg[2], const int amountOfCellsArg[2], const MeshSpacing meshSpacingArg[2], const EInitialisationMethod initialisationMethod) :
	name(name),
	initialisationMethod(initialisationMethod),
	position(position)
{
	// Not a very pretty way to do this, but initialiser lists appear to break when using c style arrays
	size[0] = sizeArg[0];
	size[1] = sizeArg[1];
	amountOfCells[0] = amountOfCellsArg[0];
	amountOfCells[1] = amountOfCellsArg[1];

	// Initialising the field quantities with the given dimension
	rho = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	u = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	v = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	p = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	E = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	T = FieldQuantity(amountOfCells[0], amountOfCells[1]);
	H = FieldQuantity(amountOfCells[0], amountOfCells[1]);

	meshSpacing[0] = MeshSpacing(meshSpacingArg[0]);
	meshSpacing[1] = MeshSpacing(meshSpacingArg[1]);

	PopulateDomainDimensions();

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

int Domain::GetTotalAmountOfCells() const
{
	return amountOfCells[0] * amountOfCells[1];
}

void Domain::GetCellSizes(const CellIndex cellPos, double& xSizeOut, double& ySizeOut) const
{
	xSizeOut = meshSpacing[0].GetCellWidth(cellPos.x);
	ySizeOut = meshSpacing[1].GetCellWidth(cellPos.y);
}

void Domain::CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to)
{
	rho.CopyToBuffer(from, to);
	u.CopyToBuffer(from, to);
	v.CopyToBuffer(from, to);
	p.CopyToBuffer(from, to);
	E.CopyToBuffer(from, to);
	T.CopyToBuffer(from, to);
	H.CopyToBuffer(from, to);
}

void Domain::SetToAmbientConditions(const double TSet, const double pSet, const double uSet, const double vSet, const double R_ideal, const double gamma)
{
	T.SetAllToValue(TSet);
	p.SetAllToValue(pSet);
	u.SetAllToValue(uSet);
	v.SetAllToValue(vSet);

	const double rhoSet = pSet / (TSet * R_ideal);
	const double ESet = pSet / (gamma - 1.0) + 0.5 * rhoSet * (pow(uSet, 2) + pow(vSet, 2));
	const double HSet = (ESet + pSet) / rhoSet;
	rho.SetAllToValue(rhoSet);
	E.SetAllToValue(ESet);
	H.SetAllToValue(HSet);
}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y!");
	}
}

void Domain::PopulateDomainDimensions()
{
	// Apologies if this function is a little hard to wrap your head around, but this saves a lot of performance pain and repeated code.

	std::vector<double> lengths[2]({ {size[0], 0}, {size[1], 0} });
	std::vector<double> centerPositions[2]({ {size[0], 0}, {size[1], 0} });
	// do it for both axes, which are iterates as meshSpacing[axis]
	for (int axis = 0; axis < 2; axis++)
	{
		// As the spacing in the y-direction is not dependent on the x position and vice versa, they can be pre-calculated, and later just populated.
		double previousPosition = 0;
		for (int i = 0; i < size[axis]; i++)
		{
			// Get the length of the current cell
			double length;
			switch (meshSpacing[axis].spacingType)
			{
			case EMeshSpacingType::CONSTANT:
				length = length / size[axis];
				break;
			case EMeshSpacingType::LINEAR:
				break;

			default:
				throw std::logic_error("Populating domain dimensions is not implemented for this mesh spacing type.");
			}
			

			lengths[axis][i] = length;
			// Note that centerPositions[i-1] cannot be used here, because this stores the center positions. On top of it, it would not be defined for the first iteration.
			centerPositions[axis][i] = previousPosition + 0.5 * length;
			previousPosition += length;
		}
	}

	// Set the x-spacing for all y cells with this x coordinate (and vice versa) as they are independent!
	//int otherAxis = (axis == 0) ? 1 : 0;
	for (int xIdx = 0; xIdx < size[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < size[1]; yIdx++)
		{
			cellLength[0].main[cellLength[0].At(xIdx, yIdx)] = lengths[0][xIdx];
			cellLength[1].main[cellLength[1].At(xIdx, yIdx)] = lengths[1][yIdx];
			localCellCenterPosition[0].main[localCellCenterPosition[0].At(xIdx, yIdx)] = centerPositions[0][xIdx];
			localCellCenterPosition[1].main[localCellCenterPosition[1].At(xIdx, yIdx)] = centerPositions[1][yIdx];
		}
	}
}

