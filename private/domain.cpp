#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cassert>


Domain::Domain(const std::string& name, SimCase* simCase, Position position, const double sizeArg[2], const int amountOfCellsArg[2], const MeshSpacing meshSpacingArg[2], const EInitialisationMethod initialisationMethod) :
	name(name),
	simCase(simCase),
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

std::pair<EBoundaryLocation, double> Domain::GetLocationAlongBoundaryInAdjacentDomain(
	const EBoundaryLocation boundaryInThisDomain, const double positionAlongBoundaryInThisDomain) const
{
	/* As it's the complement boundary, the orientation of the coordinate frame is switched too. Since the lengths of the domains and hence their boundaries are exactly the same, we can just subtract.
	 *							
	 *	Domain B				  
	 *	L	   (L-X) 			   0
	 *	+---<----+---<----<----<---+
	 *
	 *	+--->----+--->---->---->---+
	 *	0	     X				   L
	 *	
	 */

	EBoundaryLocation complement = Opposite(boundaryInThisDomain);
	double totalLength = (boundaryInThisDomain == TOP || boundaryInThisDomain== BOTTOM) ? size[1] : size[0];
	double posInOther = totalLength - positionAlongBoundaryInThisDomain;

#ifdef _DEBUG
	if (posInOther < 0 )
		throw std::logic_error("Getting the complement of this location on another boundary returned a negative value, which is invalid!");
#endif

	return std::make_pair(complement, posInOther);
}

void Domain::SetBoundaryType(const EBoundaryLocation location, const EBoundaryCondition type)
{
	if (type == EBoundaryCondition::CONNECTED)
	{
		throw std::invalid_argument("Boundary type cannot be manually set to connected. Use ConnectBoundaries() instead.");
	}
	// Use the ordering of the enum to index the boundaries[] array
	boundaries[location] = Boundary(type);
}

void Domain::ConnectBoundary(const EBoundaryLocation location, Domain* otherDomain)
{
	Boundary& thisBoundary = boundaries.at(location);
	// Validate the boundary on this domain 
	if (thisBoundary.boundaryType != NOT_SET)
		throw std::logic_error("This boundary has already been assigned another value");
	if (thisBoundary.connectedBoundary != nullptr)
		throw std::logic_error("This boundary has already been connected with another one!");

	// Get the boundary on the other domain.
	// Since we're only working with rectangular domains, a boundary will always connect to its opposite.
	Boundary& otherBoundary = otherDomain->boundaries.at(Opposite(location));

	// Validate the boundary on the other domain
	if (otherBoundary.boundaryType != NOT_SET)
		throw std::logic_error("This boundary has already been assigned another value");
	if (otherBoundary.connectedBoundary != nullptr)
		throw std::logic_error("This boundary has already been connected with another one!");

	// Make sure they're the same size and have the same amount of cells. Although the size is not a numerical necessity, it's still a physical one.
	const bool bVertical = (location == LEFT || location == RIGHT);
	if (amountOfCells[bVertical] != otherDomain->amountOfCells[bVertical])
		throw std::logic_error("Cannot connect two boundaries that have a different number of cells");
	if (!IsCloseToZero(size[bVertical] - otherDomain->size[bVertical]))
		throw std::logic_error("Cannot connect two boundaries that have a different physical size. Though numerically admissable, physically impossible.");
	if (!(meshSpacing[bVertical] == otherDomain->meshSpacing[bVertical]))
		throw std::logic_error("Cannot connect two boundaries that have a different spacing.");
	
	// All good, connect both of them.
	thisBoundary.boundaryType = CONNECTED;
	thisBoundary.connectedBoundary = &otherBoundary;
	otherBoundary.boundaryType = CONNECTED;
	otherBoundary.connectedBoundary = &thisBoundary;
}

int Domain::GetTotalAmountOfCells() const
{
	return amountOfCells[0] * amountOfCells[1];
}

CellIndex Domain::GetOriginIndexOfBoundary(const EBoundaryLocation boundary) const
{
	switch (boundary)
	{
	case BOTTOM:
		return {amountOfCells[0],0, TOP};
	case LEFT:
		return {0,0, RIGHT};
	case RIGHT:
		return {amountOfCells[0], amountOfCells[1], LEFT};
	case TOP:
		return {0, amountOfCells[1], BOTTOM};
	default:
		throw std::logic_error("GetGhostOrigin is not implemented for this boundary location.");
	}
}

std::pair<double, double> Domain::GetCellSizes(const CellIndex cellPos) const
{
	double xSizeOut = meshSpacing[0].GetCellWidth(cellPos.x);
	double ySizeOut = meshSpacing[1].GetCellWidth(cellPos.y);
	return std::make_pair(xSizeOut, ySizeOut);
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

double Domain::SpecificHeatRatio() const
{
	// todo: Extend this, so that it actually calculated based on the species and temperature in the domain.
	return 1.4;
}

void Domain::UpdateGhostCells()
{
	for (const auto& it : boundaries)
	{
		const EBoundaryLocation location = it.first();
		const Boundary& boundary = it.second;

		// Note that the below functions all work in relative coordinate frames.
		// Technically a performance gain could be achieved by not calling a transformation on the frame every time, instead directly accessing those entries directly by looping. However, this is more legible.
		switch (boundary.boundaryType) {
		case NOT_SET:
				throw std::logic_error("Cannot update Ghost cell if boundary condition is not set.");
			case SLIP:
				PopulateSlipConditionGhostCells(location);
				break;
			case NO_SLIP:
				PopulateNoSlipConditionGhostCells(location);
				break;
			case CONNECTED:
				PopulateConnectedGhostCells(location);
				break;
		default:
				throw std::logic_error("Updating ghost cells is not implemented for this type of boundary.");
		}
	}
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
	// First put it in just arrays, since the different axis are independent of eachother.
	std::vector<double> lengths[2] = { {size[0], 0}, {size[1], 0} };
	std::vector<double> centerPositions[2] = { {size[0], 0}, {size[1], 0} };
	// do it for both axes, which are iterates as meshSpacing[axis]
	for (int axis = 0; axis < 2; axis++)
	{
		// As the spacing in the y-direction is not dependent on the x position and vice versa, they can be pre-calculated, and later just populated.
		double previousPosition = 0;
		for (int i = 0; i < size[axis]; i++)
		{
			// Get the length of the current cell
			double cellLength = 0;
			switch (meshSpacing[axis].spacingType)
			{
			case EMeshSpacingType::CONSTANT:
				cellLength = cellLength / size[axis];
				break;
			case EMeshSpacingType::LINEAR:
				break;

			default:
				throw std::logic_error("Populating domain dimensions is not implemented for this mesh spacing type.");
			}
			

			lengths[axis][i] = cellLength;
			// Note that centerPositions[i-1] cannot be used here, because this stores the center positions. On top of it, it would not be defined for the first iteration.
			centerPositions[axis][i] = previousPosition + 0.5 * cellLength;
			previousPosition += cellLength;
		}
	}

	// Set the x-spacing for all y cells with this x coordinate (and vice versa) as they are independent!
	//int otherAxis = (axis == 0) ? 1 : 0;
	for (int xIdx = 0; xIdx < size[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < size[1]; yIdx++)
		{
			cellLengths[0].main[cellLengths[0].At(xIdx, yIdx)] = lengths[0][xIdx];
			cellLengths[1].main[cellLengths[1].At(xIdx, yIdx)] = lengths[1][yIdx];
			localCellCenterPositions[0].main[localCellCenterPositions[0].At(xIdx, yIdx)] = centerPositions[0][xIdx];
			localCellCenterPositions[1].main[localCellCenterPositions[1].At(xIdx, yIdx)] = centerPositions[1][yIdx];
		}
	}
}

void Domain::PopulateSlipConditionGhostCells(const EBoundaryLocation boundary)
{
	CellIndex ghostOrigin = GetOriginIndexOfBoundary(boundary);
	auto extent = GetGhostDimensions(boundary);

	// Generally constant in the local x-direction, so y outer loop.
	for (int yLocalIdx = 0; yLocalIdx < extent.second; yLocalIdx++)
	{
		for (int xLocalIdx = 0; xLocalIdx < extent.first; xLocalIdx++)
		{
			// Note that for the ghost cells, the relative location in the reference frame relative to the boundary is negative, and needs to be offset by -1 as well, as 0 is the first positive cell, and not the actual zero-line.
			const int ghostY = -yLocalIdx-1;
			const CellIndex ghostInBoundaryLocalReferenceFrame = {xLocalIdx, ghostY, boundary};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostInBoundaryLocalReferenceFrame, ghostOrigin, {0,0, TOP} );
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem({xLocalIdx, yLocalIdx, boundary}, ghostOrigin, {0,0, TOP});
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(ValidateCellIndex(sourceIndex, false));
#endif
			
			rho(ghostIndex)		=	rho(sourceIndex);
			p(ghostIndex)		=	p(sourceIndex);
			u(ghostIndex)		= - u(sourceIndex); // Flipped!
			v(ghostIndex)		=   v(sourceIndex);
			H(ghostIndex)		=	H(sourceIndex);
			E(ghostIndex)		=	E(sourceIndex);
			T(ghostIndex)		=	T(sourceIndex);
		}
	}
}

void Domain::PopulateNoSlipConditionGhostCells(const EBoundaryLocation boundary)
{
	// This almost entirely a copy-paste from Domain::PopulateSlipConditionGhostCells(), but now v is also flipped!
	CellIndex ghostOrigin = GetOriginIndexOfBoundary(boundary);
	auto extent = GetGhostDimensions(boundary);

	// Generally constant in the local x-direction, so y outer loop.
	for (int yLocalIdx = 0; yLocalIdx < extent.second; yLocalIdx++)
	{
		for (int xLocalIdx = 0; xLocalIdx < extent.first; xLocalIdx++)
		{
			// Note that for the ghost cells, the relative location in the reference frame relative to the boundary is negative, and needs to be offset by -1 as well, as 0 is the first positive cell, and not the actual zero-line.
			const int ghostY = -yLocalIdx-1;
			const CellIndex ghostInBoundaryLocalReferenceFrame = {xLocalIdx, ghostY, boundary};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostInBoundaryLocalReferenceFrame, ghostOrigin, {0,0, TOP} );
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem({xLocalIdx, yLocalIdx, boundary}, ghostOrigin, {0,0, TOP});
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(ValidateCellIndex(sourceIndex, false));
#endif

			rho(ghostIndex)		=	rho(sourceIndex);
			p(ghostIndex)		=	p(sourceIndex);
			u(ghostIndex)		= - u(sourceIndex);	// Flipped!
			v(ghostIndex)		= -	v(sourceIndex);	// Flipped!
			H(ghostIndex)		=	H(sourceIndex);
			E(ghostIndex)		=	E(sourceIndex);
			T(ghostIndex)		=	T(sourceIndex);
		}
	}
}

void Domain::PopulateConnectedGhostCells(const EBoundaryLocation boundary)
{
	/*
	 * Consider this slice for two connected domains. In this example, there are 2 ghost cells, and the domain is 1 cell thick in the axis that it is not connected to.
	 * The LEFT side of domain 1 is connected to the RIGHT side of domain 2
	 * Note that in reality, these two slices are physically in the same location, but for this drawing they have been put on top of eachother (to visualise better).
	 * Populating the ghost cells is done by taking the values from the connected domain's actual values, and putting them in the ghost cells of this domain.
	 * this means that, for domain 1: (n-2) -> G1; (n-1) -> G0
	 * conversely, for domain 2: 1 -> G1; 0 -> G2
	 * 
	 * 
	 *			Domain 2
	 * -----------------------------------------+
	 *		n-3	|	n-2	|	n-1	||	G0	|	G1	|
	 *			|		|		||		|		|
	 * ---------+-------+-------+-------+-------+-------------------
	 *			|	G1	|	G0	||	0	|	1	|	2
	 *			|		|		||		|		|
	 *			+-------+-------+-------+-------+-------------------
	 *									Domain 1
	 *
	 *					 		y	(y-datum coincides at the bottom left of cell 0).
	 *					 		|	domain reference frame
	 *					 		+ -- x
	 */

	// Validate that this is in fact a connected boundary
	if (boundaries.at(boundary).boundaryType != CONNECTED || boundaries.at(boundary).connectedBoundary == nullptr)
		throw std::logic_error("This boundary has not been connected properly.");

#ifdef _DEBUG
	// Also check the opposite boundary
	Boundary* otherBoundary = boundaries.at(boundary).connectedBoundary;
	if (otherBoundary->boundaryType != CONNECTED || otherBoundary->connectedBoundary != &boundaries.at(boundary))
		throw std::logic_error("The boundary is only properly connected one way!");
#endif

	const Domain* otherDomain = boundaries.at(boundary).connectedBoundary->domain;
	const bool bVertical = (boundary == LEFT || boundary == RIGHT);

#ifdef _DEBUG
	// Check if they have the same amount of cells.
	assert(amountOfCells[bVertical] != otherDomain->amountOfCells[bVertical]);
#endif

	/* Note that the coordinate systems' origins and positive axis directions are in opposite directions.
	 * This is both for the x coordinate and the y coordinate, as both are right-handed
	 * 
	 *	Domain 2		Domain 1
	 *
	 *  --- n - +		+ - 0 ---
	 *			|		|
	 *		2 	+		+   1
	 *			|		|
	 *		1	+		+   2
	 *			|		|
	 *	---	0 -	+		+ - n ---
	 */

	CellIndex ghostOrigin = GetOriginIndexOfBoundary(boundary);
	CellIndex complementOrigin = otherDomain->GetOriginIndexOfBoundary(Opposite(boundary));
	auto extent = GetGhostDimensions(boundary);

	// Generally constant in the local x-direction, so y outer loop.
	for (int yLocalIdx = 0; yLocalIdx < extent.second; yLocalIdx++)
	{
		for (int xLocalIdx = 0; xLocalIdx < extent.first; xLocalIdx++)
		{
			// Getting index for the ghost cell to write to
			const CellIndex ghostInBoundaryLocalReferenceFrame = {xLocalIdx, -yLocalIdx-1, boundary};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostInBoundaryLocalReferenceFrame, ghostOrigin, {0,0, TOP} );
	
			// Getting index for the cell to read from. See above ASCII sketch, as the origin is on the opposite side, subtract total amount of cells from current index to get the complement.
			int complementXIndex = otherDomain->amountOfCells[bVertical] - xLocalIdx;
			const CellIndex sourceIndexInOppositeBoundaryRelative = {complementXIndex, yLocalIdx, Opposite(boundary)};
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem(sourceIndexInOppositeBoundaryRelative, complementOrigin, {0,0, TOP});
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(otherDomain->ValidateCellIndex(sourceIndex, false));
#endif

			rho(ghostIndex)		=	otherDomain->rho.GetAt(sourceIndex);
			p(ghostIndex)		=	otherDomain->p.GetAt(sourceIndex);
			u(ghostIndex)		=   otherDomain->u.GetAt(sourceIndex);	
			v(ghostIndex)		=  	otherDomain->v.GetAt(sourceIndex);	
			H(ghostIndex)		=	otherDomain->H.GetAt(sourceIndex);
			E(ghostIndex)		=	otherDomain->E.GetAt(sourceIndex);
			T(ghostIndex)		=	otherDomain->T.GetAt(sourceIndex);
			
		}
	}
}

bool Domain::ValidateCellIndex(const CellIndex cellIndex, const bool bAllowGhostCells) const
{
	if (cellIndex.x < -nGhost*bAllowGhostCells)
		return false;
	if (cellIndex.x > amountOfCells[0] + nGhost*bAllowGhostCells)
		return false;
	if (cellIndex.y < -nGhost*bAllowGhostCells)
		return false;
	if (cellIndex.y > amountOfCells[1] + nGhost*bAllowGhostCells)
		return false;
	
	return true;
}

std::pair<int, int> Domain::GetGhostDimensions(EBoundaryLocation boundary)
{
	// Note that these sorta flip the coordinate system!
	switch (boundary)
	{
	case BOTTOM:
		return {amountOfCells[0],nGhost};
	case LEFT:
		return {amountOfCells[1],nGhost};
	case RIGHT:
		return {amountOfCells[1],nGhost};
	case TOP:
		return {amountOfCells[0],nGhost};
	default:
		throw std::logic_error("GetGhostOrigin is not implemented for this boundary location.");
	}
}

