#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cassert>

#include "euler_container.h"
#include "flux_splitting.h"


Domain::Domain(const std::string& name, SimCase* simCase, const Position& position, const double sizeArg[2], const int amountOfCellsArg[2], const MeshSpacing meshSpacingArg[2], const EInitialisationMethod initialisationMethod, const int ghostCellDepth) :
	name(name),
	simCase(simCase),
	initialisationMethod(initialisationMethod),
	position(position),
	nGhost(ghostCellDepth)
{
	// Not a very pretty way to do this, but initialiser lists appear to break when using c style arrays
	size[0] = sizeArg[0];
	size[1] = sizeArg[1];
	amountOfCells[0] = amountOfCellsArg[0];
	amountOfCells[1] = amountOfCellsArg[1];

	// Initialising the field quantities with the given dimension
	rho = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	u = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	v = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	p = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	E = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	T = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	H = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);

	meshSpacing[0] = MeshSpacing(meshSpacingArg[0]);
	meshSpacing[1] = MeshSpacing(meshSpacingArg[1]);

	CacheCellSizes();

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

void Domain::SetToAmbientConditions(const double temperatureSet, const double pSet, const double uSet, const double vSet, const double R_ideal, const double gamma)
{
	T.currentTimeStep.SetAllToValue(temperatureSet);
	p.currentTimeStep.SetAllToValue(pSet);
	u.currentTimeStep.SetAllToValue(uSet);
	v.currentTimeStep.SetAllToValue(vSet);

	const double rhoSet = pSet / (temperatureSet * R_ideal);
	const double ESet = pSet / (gamma - 1.0) + 0.5 * rhoSet * (pow(uSet, 2) + pow(vSet, 2));
	const double HSet = (ESet + pSet) / rhoSet;
	rho.currentTimeStep.SetAllToValue(rhoSet);
	E.currentTimeStep.SetAllToValue(ESet);
	H.currentTimeStep.SetAllToValue(HSet);
}

double Domain::SpecificHeatRatio() const
{
	// todo: Extend this, so that it actually calculated based on the species and temperature in the domain.
	return 1.4;
}

double Domain::GasConstant() const
{
	return 287.0;
}

void Domain::UpdateGhostCells()
{
	// These can all be parallelised!
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

void Domain::PopulateFlowDeltaBuffer(const double dt, const double AUSMkFactor, const double entropyFix)
{
	// Depending on whether or not there are shock fronts, the integration scheme changes. Determining whether or not there are shock fronts is done by looking at the MUSCL interpolated values of the field quantities. So, first calculate & cache the MUSCL vales.
	
	const double bias = simCase->MUSCLBias;
	const EFluxLimiterType fluxLimiter = simCase->fluxLimiterType;
	const int rungeKuttaOrder = simCase->rungeKuttaOrder;

	rho.PopulateMUSCLBuffers(RUNGE_KUTTA, bias, fluxLimiter);
	u.PopulateMUSCLBuffers(RUNGE_KUTTA, bias, fluxLimiter);
	v.PopulateMUSCLBuffers(RUNGE_KUTTA, bias, fluxLimiter);
	p.PopulateMUSCLBuffers(RUNGE_KUTTA, bias, fluxLimiter);

	// Do flux splitting on all the faces. Use hanel if flow is sonic in either cell, ausm_dv if not.
	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			// These here contain 4 terms for the euler equations. They all represent the flux of those variables at each side of the cell.
			EulerContinuity leftFlux, rightFlux, upFlux, downFlux;
			
			// rho, u, v and p's MUSCL values are stored in their respective buffers. the other variables are not, but they can be calculated directly from the state values in those buffers. So, do this for energy and enthalpy.
			double gamma = SpecificHeatRatio();
			double gasConstant = GasConstant();
			// energy
			const double eLeft = p.leftFaceMUSCLBuffer.GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.leftFaceMUSCLBuffer.GetAt(xIdx, yIdx) * (pow(u.leftFaceMUSCLBuffer(xIdx, yIdx), 2) + pow(v.leftFaceMUSCLBuffer(xIdx, yIdx), 2));
			const double eRight = p.rightFaceMUSCLBuffer.GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.rightFaceMUSCLBuffer.GetAt(xIdx, yIdx) * (pow(u.rightFaceMUSCLBuffer(xIdx, yIdx), 2) + pow(v.rightFaceMUSCLBuffer(xIdx, yIdx), 2));
			const double eTop = p.topFaceMUSCLBuffer.GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.topFaceMUSCLBuffer.GetAt(xIdx, yIdx) * (pow(u.topFaceMUSCLBuffer(xIdx, yIdx), 2) + pow(v.topFaceMUSCLBuffer(xIdx, yIdx), 2));
			const double eBottom = p.bottomFaceMUSCLBuffer.GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.bottomFaceMUSCLBuffer.GetAt(xIdx, yIdx) * (pow(u.bottomFaceMUSCLBuffer(xIdx, yIdx), 2) + pow(v.bottomFaceMUSCLBuffer(xIdx, yIdx), 2));

			// enthalpy
			const double hLeft = (eLeft + p.leftFaceMUSCLBuffer.GetAt(xIdx,yIdx))/rho.leftFaceMUSCLBuffer.GetAt(xIdx,yIdx);
			const double hRight = (eRight + p.rightFaceMUSCLBuffer.GetAt(xIdx,yIdx))/rho.rightFaceMUSCLBuffer.GetAt(xIdx,yIdx);
			const double hTop = (eTop + p.topFaceMUSCLBuffer.GetAt(xIdx,yIdx))/rho.topFaceMUSCLBuffer.GetAt(xIdx,yIdx);
			const double hBottom = (eBottom + p.bottomFaceMUSCLBuffer.GetAt(xIdx,yIdx))/rho.bottomFaceMUSCLBuffer.GetAt(xIdx,yIdx);

			// Cache whether or not the flow is supersonic through faces.
			const double cLeft = sqrt(gamma * p.leftFaceMUSCLBuffer.GetAt(xIdx, yIdx)/rho.leftFaceMUSCLBuffer.GetAt(xIdx, yIdx));
			const double cRight = sqrt(gamma * p.rightFaceMUSCLBuffer.GetAt(xIdx, yIdx)/rho.rightFaceMUSCLBuffer.GetAt(xIdx, yIdx));
			const double cTop = sqrt(gamma * p.topFaceMUSCLBuffer.GetAt(xIdx, yIdx)/rho.topFaceMUSCLBuffer.GetAt(xIdx, yIdx));
			const double cBottom = sqrt(gamma * p.bottomFaceMUSCLBuffer.GetAt(xIdx, yIdx)/rho.bottomFaceMUSCLBuffer.GetAt(xIdx, yIdx));
			const double speedOfSound[4] = {cRight, cLeft, cTop, cBottom};

			// Now, for every face, determine which flux scheme to use based on whether or not it is critical (sonic), and perform the flux transfer.
			// indexes for sides ---- 0: right; 1: left; 2: top; 3: down;
			EBoundaryLocation faces[4] = {LEFT, RIGHT, TOP, BOTTOM};
			EulerContinuity fluxSplit[4];
			for (int sideIndex = 0; sideIndex < 4; sideIndex++)
			{
				EBoundaryLocation face = faces[sideIndex]; // Easier to debug and read, just a map of sideIndex. Note that enum decays into int.
				CellIndex rf; // The position in the MUSCL buffer where the relevant values are stored.
				// the fluxes for the cells that are 'backwards' in the positive axis are therefore by definition equal to the

				/* The values in the MUSCL buffers are stored as positive terms on their current cells; this means that musclbuffer(i,j) is the flux to the positive axis.
				 * There are 4 buffers; each for flow going in one specific direction.
				 * Hence, for flow in the right direction, take the flux in the right direction at right face, and the left face; right face is stored at i,j normally, and for the one on the left face use the value for the previous cell (identical).
				 *	
				 *									|	cell (xIdx, yIdx)	|
				 *		musclRight(xIdx - 1, yIdx) -|-->				  --|--> musclRight(xIdx, yIdx)
				 *									|						|
				 *
				 */
				
				if (face == LEFT)
					rf = CellIndex(xIdx -1, yIdx);
				if (sideIndex == BOTTOM)
					rf = CellIndex(xIdx, yIdx - 1);

				EulerContinuity leftContinuity, rightContinuity;
				if (face == LEFT || face == RIGHT) // horizontal flux
				{
					leftContinuity = EulerContinuity(rho.leftFaceMUSCLBuffer.GetAt(rf), u.leftFaceMUSCLBuffer.GetAt(rf), v.leftFaceMUSCLBuffer.GetAt(rf), eLeft, hLeft, p.leftFaceMUSCLBuffer.GetAt(rf));
					rightContinuity = EulerContinuity(rho.rightFaceMUSCLBuffer.GetAt(rf), u.rightFaceMUSCLBuffer.GetAt(rf), v.rightFaceMUSCLBuffer.GetAt(rf), eRight, hRight, p.rightFaceMUSCLBuffer.GetAt(rf));
				}
				else //face == TOP || face == BOTTOM, vertical flux
				{
					leftContinuity = EulerContinuity(rho.leftFaceMUSCLBuffer.GetAt(rf), u.leftFaceMUSCLBuffer.GetAt(rf), v.leftFaceMUSCLBuffer.GetAt(rf), E.leftFaceMUSCLBuffer.GetAt(rf), H.leftFaceMUSCLBuffer.GetAt(rf), p.leftFaceMUSCLBuffer.GetAt(rf));
					rightContinuity = EulerContinuity(rho.rightFaceMUSCLBuffer.GetAt(rf), u.rightFaceMUSCLBuffer.GetAt(rf), v.rightFaceMUSCLBuffer.GetAt(rf), E.rightFaceMUSCLBuffer.GetAt(rf), H.rightFaceMUSCLBuffer.GetAt(rf), p.rightFaceMUSCLBuffer.GetAt(rf));
				}

				// Set the flux split for this side (declared above in the array).
				if (v.rightFaceMUSCLBuffer.GetAt(rf) > speedOfSound[sideIndex])
				{
					// It's sonic, use Hanel.
					fluxSplit[sideIndex] = HanelFluxSplitting(leftContinuity, rightContinuity, gamma, AUSMkFactor, entropyFix);
				}
				else
				{
					// It's subsonic, use AUSM_DV.
					fluxSplit[sideIndex] = AUSMDVFluxSplitting(leftContinuity, rightContinuity, gasConstant, gamma, AUSMkFactor, entropyFix);
				}
			}

			// Total accumulation is what goes in - what goes out
			const CellIndex currentCell(xIdx, yIdx);
			auto cellSizes = GetCellSizes(currentCell);
			//EulerContinuity accumulation = ((rightFlux - leftFlux)/cellSizes.first + (upFlux - downFlux)/cellSizes.second)*dt; // dy = dr in this case, if you compensate for the squashification.
			rho.deltaDueToFlow(xIdx,yIdx) = ((fluxSplit[RIGHT].density - fluxSplit[LEFT].density)/cellSizes.first + (fluxSplit[TOP].density - fluxSplit[BOTTOM].density)/cellSizes.second)*dt;
			v.deltaDueToFlow(xIdx,yIdx) = ((fluxSplit[RIGHT].v - fluxSplit[LEFT].v)/cellSizes.first + (fluxSplit[TOP].v - fluxSplit[BOTTOM].v)/cellSizes.second)*dt;
			u.deltaDueToFlow(xIdx,yIdx) = ((fluxSplit[RIGHT].u - fluxSplit[LEFT].u)/cellSizes.first + (fluxSplit[TOP].u - fluxSplit[BOTTOM].u)/cellSizes.second)*dt;
			E.deltaDueToFlow(xIdx,yIdx) = ((fluxSplit[RIGHT].e - fluxSplit[LEFT].e)/cellSizes.first + (fluxSplit[TOP].e - fluxSplit[BOTTOM].e)/cellSizes.second)*dt;

			// Extra term to compensate for the fact that the cell sizes are not uniform because we are in a cylindrical coordinate system. In Florian (2017), this is the term *** Hr ***.
			EulerContinuity hr;
			const double density = rho.currentTimeStep.GetAt(currentCell);
			const double xVel = u.currentTimeStep.GetAt(currentCell);
			const double yVel = v.currentTimeStep.GetAt(currentCell);
			const double yc = localCellCenterPositions[1].GetAt(currentCell);
			const double enthalpy = H.currentTimeStep.GetAt(currentCell);
			
			rho.deltaDueToFlow(xIdx,yIdx) += density * yVel / yc * dt;
			u.deltaDueToFlow(xIdx,yIdx) = density * yVel * xVel / yc * dt;
			v.deltaDueToFlow(xIdx,yIdx) =  density * yVel * yVel / yc * dt;
			E.deltaDueToFlow(xIdx,yIdx) =  density * yVel * enthalpy / yc * dt;
		}
	}	
}

void Domain::SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(const int currentRungeKuttaIter)
{
#ifdef _DEBUG
	// Fail if the buffers are empty; they clearly must be set.
	assert(!rho.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(!u.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(!v.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(!E.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(rho.deltaDueToFlow.IsFilledWithZeroes());
	assert(rho.deltaDueToValve.IsFilledWithZeroes());
	assert(u.deltaDueToFlow.IsFilledWithZeroes());
	assert(u.deltaDueToValve.IsFilledWithZeroes());
	assert(v.deltaDueToFlow.IsFilledWithZeroes());
	assert(v.deltaDueToValve.IsFilledWithZeroes());
	assert(E.deltaDueToFlow.IsFilledWithZeroes());
	assert(E.deltaDueToValve.IsFilledWithZeroes());
#endif
	
	const int rungeKuttaOrder = simCase->rungeKuttaOrder;
	double rkK = 1./(rungeKuttaOrder-currentRungeKuttaIter); // Runge-kutta factor

	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			// These are state variables, and are explicitly expressed. combine the values in the delta buffers, and apply runge-kutta scaling.
			double dRho = rho.deltaDueToFlow.GetAt(xIdx, yIdx) + rho.deltaDueToValve.GetAt(xIdx, yIdx);
			double dU = u.deltaDueToFlow.GetAt(xIdx, yIdx) + u.deltaDueToValve.GetAt(xIdx, yIdx);
			double dV = v.deltaDueToFlow.GetAt(xIdx, yIdx) + v.deltaDueToValve.GetAt(xIdx, yIdx);
			double dE = E.deltaDueToFlow.GetAt(xIdx, yIdx) + E.deltaDueToValve.GetAt(xIdx, yIdx);
			rho.nextTimeStepBuffer(xIdx, yIdx) = rho.currentTimeStep.GetAt(xIdx, yIdx) - rkK * dRho;
			u.nextTimeStepBuffer(xIdx, yIdx) = (rho.currentTimeStep.GetAt(xIdx, yIdx) * u.currentTimeStep.GetAt(xIdx,yIdx) - rkK * dU) / rho.nextTimeStepBuffer(xIdx, yIdx);
			v.nextTimeStepBuffer(xIdx, yIdx) = (rho.currentTimeStep.GetAt(xIdx, yIdx) * v.currentTimeStep.GetAt(xIdx,yIdx) - rkK * dV) / rho.nextTimeStepBuffer(xIdx, yIdx);
			E.nextTimeStepBuffer(xIdx, yIdx) = E.currentTimeStep.GetAt(xIdx, yIdx) - rkK * dRho;

			// The others are not state variables; they can be calculated using the known variables. Calculate them now.
			//todo: set a build mode where these are not calculated to speed things up, as this is only really necessary for data export.
			p.nextTimeStepBuffer(xIdx, yIdx) = (SpecificHeatRatio()-1) * E.nextTimeStepBuffer.GetAt(xIdx, yIdx) - 0.5*rho.nextTimeStepBuffer(xIdx,yIdx)*pow(u.nextTimeStepBuffer.GetAt(xIdx,yIdx) + v.nextTimeStepBuffer.GetAt(xIdx,yIdx), 2);
			T.nextTimeStepBuffer(xIdx, yIdx) = p.nextTimeStepBuffer.GetAt(xIdx, yIdx) / (GasConstant() * rho.nextTimeStepBuffer.GetAt(xIdx,yIdx));
			H.nextTimeStepBuffer(xIdx,yIdx) = (E.nextTimeStepBuffer.GetAt(xIdx,yIdx) + p.nextTimeStepBuffer(xIdx,yIdx))/rho.nextTimeStepBuffer.GetAt(xIdx,yIdx);
		}
	}

	// Clean up, empty the buffers
	rho.deltaDueToFlow.SetAllToValue(0);
	rho.deltaDueToValve.SetAllToValue(0);
	u.deltaDueToFlow.SetAllToValue(0);
	u.deltaDueToValve.SetAllToValue(0);
	v.deltaDueToFlow.SetAllToValue(0);
	v.deltaDueToValve.SetAllToValue(0);
	E.deltaDueToFlow.SetAllToValue(0);
	E.deltaDueToValve.SetAllToValue(0);
	// Note that clearing the others (non-state variables) is not necessary, as they are not written to.
}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y!");
	}
}

void Domain::CacheCellSizes()
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
			cellLengths[0](xIdx,yIdx) = lengths[0].at(xIdx);
			cellLengths[1](xIdx,yIdx) = lengths[1].at(yIdx);
			localCellCenterPositions[0](xIdx, yIdx) = centerPositions[0].at(xIdx);
			localCellCenterPositions[1](xIdx, yIdx) = centerPositions[1].at(yIdx);
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
			
			rho.rungeKuttaBuffer(ghostIndex)		=	rho.rungeKuttaBuffer(sourceIndex);
			p.rungeKuttaBuffer(ghostIndex)		=	p.rungeKuttaBuffer(sourceIndex);
			u.rungeKuttaBuffer(ghostIndex)		= - u.rungeKuttaBuffer(sourceIndex); // Flipped!
			v.rungeKuttaBuffer(ghostIndex)		=   v.rungeKuttaBuffer(sourceIndex);
			H.rungeKuttaBuffer(ghostIndex)		=	H.rungeKuttaBuffer(sourceIndex);
			E.rungeKuttaBuffer(ghostIndex)		=	E.rungeKuttaBuffer(sourceIndex);
			T.rungeKuttaBuffer(ghostIndex)		=	T.rungeKuttaBuffer(sourceIndex);
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

			rho.rungeKuttaBuffer(ghostIndex)	=	rho.rungeKuttaBuffer(sourceIndex);
			p.rungeKuttaBuffer(ghostIndex)		=	p.rungeKuttaBuffer(sourceIndex);
			u.rungeKuttaBuffer(ghostIndex)		= - u.rungeKuttaBuffer(sourceIndex);	// Flipped!
			v.rungeKuttaBuffer(ghostIndex)		= -	v.rungeKuttaBuffer(sourceIndex);	// Flipped!
			H.rungeKuttaBuffer(ghostIndex)		=	H.rungeKuttaBuffer(sourceIndex);
			E.rungeKuttaBuffer(ghostIndex)		=	E.rungeKuttaBuffer(sourceIndex);
			T.rungeKuttaBuffer(ghostIndex)		=	T.rungeKuttaBuffer(sourceIndex);
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

			rho.rungeKuttaBuffer(ghostIndex)		=	otherDomain->rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer(ghostIndex)		=	otherDomain->p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer(ghostIndex)		=   otherDomain->u.rungeKuttaBuffer.GetAt(sourceIndex);	
			v.rungeKuttaBuffer(ghostIndex)		=  	otherDomain->v.rungeKuttaBuffer.GetAt(sourceIndex);	
			H.rungeKuttaBuffer(ghostIndex)		=	otherDomain->H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer(ghostIndex)		=	otherDomain->E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer(ghostIndex)		=	otherDomain->T.rungeKuttaBuffer.GetAt(sourceIndex);
			
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

