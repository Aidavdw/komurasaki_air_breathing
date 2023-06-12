#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cassert>
#include <ios>
#include <iostream>

#include "euler_container.h"
#include "flux_splitting.h"


Domain::Domain(const std::string& name, SimCase* simCase, const Position& position, const std::pair<double, double> sizeArg, const std::pair<MeshSpacing, MeshSpacing> meshSpacingArg, const EInitialisationMethod initialisationMethod, const int ghostCellDepth) :
	name(name),
	simCase(simCase),
	initialisationMethod(initialisationMethod),
	position(position),
	nGhost(ghostCellDepth)
{
	// Not a very pretty way to do this, but initialiser lists appear to break when using c style arrays
	size[0] = sizeArg.first;
	size[1] = sizeArg.second;
	amountOfCells[0] = meshSpacingArg.first.amountOfElements;
	amountOfCells[1] = meshSpacingArg.second.amountOfElements;

	// Initialising the field quantities with the given dimension
	rho = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	u = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	v = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	p = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	E = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	T = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);
	H = FieldQuantity(this, amountOfCells[0], amountOfCells[1], 0, nGhost);

	// todo: change how this is done, because right now the length & resolution need to be input twice. Maybe make a proxy constructor?
	meshSpacing[0] = MeshSpacing(meshSpacingArg.first);
	meshSpacing[1] = MeshSpacing(meshSpacingArg.second);

	CacheCellSizes();

}

CellIndex Domain::InvertPositionToIndex(const Position pos, Position& distanceFromCenterOut) const
{
	CellIndex out;
	// Do a bisection between centre positions on both axes until the distance in steps is equal to one.
	// todo: implement are more efficient method than bisection.
	{
		const size_t leftFromIndex = FindIndexLeftOfValueByBisection(localCellCenterPositions[0], pos.x);
		if (IsCloserToLeftThanToRight(pos.x, localCellCenterPositions[0].at(leftFromIndex), localCellCenterPositions[0].at(leftFromIndex+1 )))
		{
			out.x = static_cast<int>(leftFromIndex);
			distanceFromCenterOut.x = pos.x - localCellCenterPositions[0].at(leftFromIndex);
		}
		else
		{
			out.x = static_cast<int>(leftFromIndex) + 1;
			distanceFromCenterOut.x = pos.x - localCellCenterPositions[0].at(leftFromIndex + 1);
		}
	}
	{
		const size_t leftFromIndex = FindIndexLeftOfValueByBisection(localCellCenterPositions[1], pos.y);
		if (IsCloserToLeftThanToRight(pos.y, localCellCenterPositions[0].at(leftFromIndex), localCellCenterPositions[0].at(leftFromIndex+1 )))
		{
			out.y = out.x = static_cast<int>(leftFromIndex);
			distanceFromCenterOut.y = pos.y - localCellCenterPositions[0].at(leftFromIndex);
		}
		else
		{
			out.y = static_cast<int>(leftFromIndex) + 1;
			distanceFromCenterOut.y = pos.y - localCellCenterPositions[0].at(leftFromIndex + 1);
		}
	}

	return out;
}

std::pair<EFace, double> Domain::GetLocationAlongBoundaryInAdjacentDomain(
	const EFace boundaryInThisDomain, const double positionAlongBoundaryInThisDomain) const
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

	EFace complement = Opposite(boundaryInThisDomain);
	double totalLength = (boundaryInThisDomain == TOP || boundaryInThisDomain== BOTTOM) ? size[1] : size[0];
	double posInOther = totalLength - positionAlongBoundaryInThisDomain;

#ifdef _DEBUG
	if (posInOther < 0 )
		throw std::logic_error("Getting the complement of this location on another boundary returned a negative value, which is invalid!");
#endif

	return std::make_pair(complement, posInOther);
}

Position Domain::PositionAlongBoundaryToCoordinate(const EFace boundary, const double positionAlongBoundary,
	const double depth) const
{

	/*
	*				    TOP
	*					/\
	*					|
	*		 + -- > --- > --- > --- +
	*	L	 |			 			|		R
	*	E	/\			 			/\		I
	*	F	 | ->					|  ->	G
	*	T	/\          /\			/\		H
	*		 |          |			|		T
	*		 + -- > --- > --- > --- +
	*				  BOTTOM
	*/

	switch (boundary)
	{
	case BOTTOM:
		return {positionAlongBoundary, depth};
	case TOP:
		return {positionAlongBoundary, size[1] - depth};
	case LEFT:
		return {depth, positionAlongBoundary};
	case RIGHT:
		return {size[0] - depth, positionAlongBoundary};
	default:
			throw std::logic_error("PositionAlongBoundaryToCoordinate is not implemented for this face.");
	}
}

void Domain::SetBoundaryType(const EFace location, const EBoundaryCondition type)
{
	if (type == EBoundaryCondition::CONNECTED)
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

CellIndex Domain::GetOriginIndexOfBoundary(const EFace boundary) const
{
	// Fixed for new coordinate system
	
	/*
	*				    TOP
	*					/\
	*					|
	*		 + -- > --- > --- > --- +
	*	L	 |			 			|		R
	*	E	/\			 			/\		I
	*	F	 | ->					|  ->	G
	*	T	/\          /\			/\		H
	*		 |          |			|		T
	*		 + -- > --- > --- > --- +
	*				  BOTTOM
	*/
	
	switch (boundary)
	{
	case BOTTOM:
		return {0,0, BOTTOM};
	case LEFT:
		return {0,0, BOTTOM};
	case RIGHT:
		return {amountOfCells[0] -1, 0, BOTTOM};
	case TOP:
		return {0, amountOfCells[1] -1, BOTTOM};
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

double Domain::GetCellVolume(const CellIndex cix) const
{
	double cellDX = cellLengths[0].at(cix.x);
	double cellDR = cellLengths[1].at(cix.y);
	double cellR = localCellCenterPositions[1].at(cix.y);
	double rInner = cellR - 0.5*cellDR;
	double rOuter = cellR + 0.5*cellDR;
	double volume = cellDX * M_PI * (pow(rOuter, 2) - pow(rInner, 2));
	return volume;
}

double Domain::GetLengthOfSide(const EFace face) const
{
	switch (face)
	{
	case LEFT: return size[1];
	case RIGHT: return size[1];
	case TOP: return size[0];
	case BOTTOM: return size[0];
	default:
		throw std::logic_error("Cannot get length of this type of side.");
	}
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

void Domain::SetToAmbientConditions(const double temperatureSet, const double pSet, const double uSet, const double vSet)
{
	T.currentTimeStep.SetAllToValue(temperatureSet);
	p.currentTimeStep.SetAllToValue(pSet);
	u.currentTimeStep.SetAllToValue(uSet);
	v.currentTimeStep.SetAllToValue(vSet);

	// Since it's uniform, gamma and ideal gas constant are the same everywhere
	const double gamma = SpecificHeatRatio();
	const double R = GasConstant();

	const double rhoSet = pSet / (temperatureSet * R);
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
	// These can all be parallelized!
	for (const auto& it : boundaries)
	{
		const EFace location = it.first;
		const Boundary& boundary = it.second;

#ifdef _DEBUG
		std::cout << "Updating Ghost Cells for domain '" << name << "', side" << FaceToString(location) << std::endl;
#endif

		// Note that the below functions all work in relative coordinate frames.
		// Technically a performance gain could be achieved by not calling a transformation on the frame every time, instead directly accessing those entries directly by looping. However, this is more legible.
		switch (boundary.boundaryType) {
		case EBoundaryCondition::NOT_SET:
				throw std::logic_error("Cannot update Ghost cell if boundary condition is not set.");
			case EBoundaryCondition::SLIP:
				PopulateSlipConditionGhostCells(location);
				break;
			case EBoundaryCondition::NO_SLIP:
				PopulateNoSlipConditionGhostCells(location);
				break;
			case EBoundaryCondition::CONNECTED:
				PopulateConnectedGhostCells(location);
				break;
		default:
				throw std::logic_error("Updating ghost cells is not implemented for this type of boundary.");
		}
	}
}

void Domain::PopulateFlowDeltaBuffer(const double dt)
{
	// Depending on whether or not there are shock fronts, the integration scheme changes. Determining whether or not there are shock fronts is done by looking at the MUSCL interpolated values of the field quantities. So, first calculate & cache the MUSCL vales.
	const SolverSettings& solverSettings = simCase->solverSettings;

#ifdef _DEBUG
	std::cout << "Calculating flow delta, and populating Flow Delta buffer of domain '" << name << "'" << std::endl;
	
	// Ensure that the buffers are actually empty to start with
	if (!rho.flux.IsFilledWithZeroes())
		throw std::logic_error("rho delta buffer is non-empty");
	if (!u.flux.IsFilledWithZeroes())
		throw std::logic_error("u delta buffer is non-empty");
	if (!v.flux.IsFilledWithZeroes())
		throw std::logic_error("v delta buffer is non-empty");
	if (!E.flux.IsFilledWithZeroes())
		throw std::logic_error("E delta buffer is non-empty");

	std::cout << "Calculating MUSCL interpolation values for rho, u, v, p" << std::endl;
#endif

	//todo: investigate how this works if the cells are not equally sized; how does accumulation work in the ghost cells?

	rho.PopulateMUSCLBuffers(EFieldQuantityBuffer::RUNGE_KUTTA, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
	u.PopulateMUSCLBuffers(EFieldQuantityBuffer::RUNGE_KUTTA, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
	v.PopulateMUSCLBuffers(EFieldQuantityBuffer::RUNGE_KUTTA, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
	p.PopulateMUSCLBuffers(EFieldQuantityBuffer::RUNGE_KUTTA, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);

	// Do flux splitting on all the faces. Use hanel if flow is sonic in either cell, ausm_dv if not.
	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			// rho, u, v and p's MUSCL values are stored in their respective buffers. the other variables are not, but they can be calculated directly from the state values in those buffers. So, do this for energy and enthalpy.
			double gamma = SpecificHeatRatio();
			double gasConstant = GasConstant();
			// energy
			const double eLeft = p.MUSCLBuffer[LEFT].GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.MUSCLBuffer[LEFT].GetAt(xIdx, yIdx) * (pow(u.MUSCLBuffer[LEFT](xIdx, yIdx), 2) + pow(v.MUSCLBuffer[LEFT](xIdx, yIdx), 2));
			const double eRight = p.MUSCLBuffer[RIGHT].GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.MUSCLBuffer[RIGHT].GetAt(xIdx, yIdx) * (pow(u.MUSCLBuffer[RIGHT](xIdx, yIdx), 2) + pow(v.MUSCLBuffer[RIGHT](xIdx, yIdx), 2));
			const double eTop = p.MUSCLBuffer[TOP].GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.MUSCLBuffer[TOP].GetAt(xIdx, yIdx) * (pow(u.MUSCLBuffer[TOP](xIdx, yIdx), 2) + pow(v.MUSCLBuffer[TOP](xIdx, yIdx), 2));
			const double eBottom = p.MUSCLBuffer[BOTTOM].GetAt(xIdx,yIdx) / (gamma - 1) + 0.5 * rho.MUSCLBuffer[BOTTOM].GetAt(xIdx, yIdx) * (pow(u.MUSCLBuffer[BOTTOM](xIdx, yIdx), 2) + pow(v.MUSCLBuffer[BOTTOM](xIdx, yIdx), 2));

			// enthalpy
			const double hLeft = (eLeft + p.MUSCLBuffer[LEFT].GetAt(xIdx,yIdx))/rho.MUSCLBuffer[LEFT].GetAt(xIdx,yIdx);
			const double hRight = (eRight + p.MUSCLBuffer[RIGHT].GetAt(xIdx,yIdx))/rho.MUSCLBuffer[RIGHT].GetAt(xIdx,yIdx);
			const double hTop = (eTop + p.MUSCLBuffer[TOP].GetAt(xIdx,yIdx))/rho.MUSCLBuffer[TOP].GetAt(xIdx,yIdx);
			const double hBottom = (eBottom + p.MUSCLBuffer[BOTTOM].GetAt(xIdx,yIdx))/rho.MUSCLBuffer[BOTTOM].GetAt(xIdx,yIdx);

			// Now, for every face, determine which flux scheme to use based on whether or not it is critical (sonic), and perform the flux transfer.
			EFace faces[4] = {LEFT, RIGHT, TOP, BOTTOM}; // indexes for sides ---- 0: right; 1: left; 2: top; 3: down;
			EulerContinuity fluxSplit[4]; // These here contain 4 terms for the euler equations. They all represent the flux of those variables at each side of the cell.
			for (int sideIndex = 0; sideIndex < 4; sideIndex++)
			{
				EFace face = faces[sideIndex]; // Easier to debug and read, just a map of sideIndex. Note that enum decays into int.
				CellIndex rf; // The position in the MUSCL buffer where the relevant values are stored.
				// the fluxes for the cells that are 'backwards' in the positive axis are therefore by definition equal to the

				/* The values in the MUSCL buffers are stored as positive terms on their current cells; this means that musclbuffer(i,j) is the flux to the positive axis.
				 * There are 4 buffers; each for flow going in one specific direction.
				 * Hence, for flow in the right direction, take the flux in the right direction at right face, and the left face; right face is stored at i,j normally, and for the one on the left face use the value for the previous cell (identical).
				 *	
				 *									|	cell (xIdx, yIdx)	|
				 *		musclRight(xIdx - 1, yIdx) -|-->				  --|--> musclRight(xIdx, yIdx) (positive)
				 *									|						|
				 *		musclLeft(xIdx -1, yIdx) <--|--					 <--|-- musclLeft(xIdx, yIdx) (negative)
				 *									|						|
				 *
				 */

				const double speedOfSound = sqrt(gamma * p.MUSCLBuffer[face].GetAt(xIdx, yIdx)/rho.MUSCLBuffer[face].GetAt(xIdx, yIdx));

				// I Don't remember what I did here, I think this can be ignored?
				/*
				if (face == LEFT)
					rf = CellIndex(xIdx -1, yIdx);
				if (face == BOTTOM)
					rf = CellIndex(xIdx, yIdx - 1);
				*/

				// Get the left- and right bound fluxes at the specific face we're considering now.
				EulerContinuity continuityInNegativeDirection;
				const EFace anti = Opposite(face);
				continuityInNegativeDirection.density = rho.MUSCLBuffer[anti].GetAt(rf);
				continuityInNegativeDirection.u = u.MUSCLBuffer[anti].GetAt(rf);
				continuityInNegativeDirection.v = v.MUSCLBuffer[anti].GetAt(rf);
				continuityInNegativeDirection.e = eLeft;
				continuityInNegativeDirection.h = hLeft;
				continuityInNegativeDirection.p = p.MUSCLBuffer[anti].GetAt(rf);
				
				EulerContinuity continuityPositiveDirection;
				continuityPositiveDirection.density = rho.MUSCLBuffer[face].GetAt(rf);
				continuityPositiveDirection.u = u.MUSCLBuffer[face].GetAt(rf);
				continuityPositiveDirection.v = v.MUSCLBuffer[face].GetAt(rf);
				continuityPositiveDirection.e = eRight;
				continuityPositiveDirection.h = hLeft;
				continuityPositiveDirection.p = p.MUSCLBuffer[face].GetAt(rf);

				// Set the flux split for this side (declared above in the array).
				if (v.MUSCLBuffer[RIGHT].GetAt(rf) > speedOfSound)
				{
					// It's sonic, use Hanel.
					fluxSplit[face] = HanelFluxSplitting(continuityInNegativeDirection, continuityPositiveDirection, gamma, solverSettings.entropyFix);
				}
				else
				{
					// It's subsonic, use AUSM_DV.
					fluxSplit[face] = AUSMDVFluxSplitting(continuityInNegativeDirection, continuityPositiveDirection, gamma, solverSettings.AUSMSwitchBias, solverSettings.entropyFix);
				}

				// If it's a vertical flux, then the u and v are in the local reference frame, and hence they must be inverted.
				if (face == TOP || face == BOTTOM)
				{
					double buf = fluxSplit[face].u;
					fluxSplit[face].u = fluxSplit[face].v;
					fluxSplit[face].v = buf;
				}
			}

			// Total accumulation is what goes in - what goes out
			const CellIndex currentCell(xIdx, yIdx);
			auto cellSizes = GetCellSizes(currentCell);
			//EulerContinuity accumulation = ((rightFlux - leftFlux)/cellSizes.first + (upFlux - downFlux)/cellSizes.second)*dt; // dy = dr in this case, if you compensate for the squashification.
			const double rhoFlux = ((fluxSplit[RIGHT].density - fluxSplit[LEFT].density)/cellSizes.first + (fluxSplit[TOP].density - fluxSplit[BOTTOM].density)/cellSizes.second)*dt;
			const double vFlux = ((fluxSplit[RIGHT].v - fluxSplit[LEFT].v)/cellSizes.first + (fluxSplit[TOP].v - fluxSplit[BOTTOM].v)/cellSizes.second)*dt;
			const double uFlux = ((fluxSplit[RIGHT].u - fluxSplit[LEFT].u)/cellSizes.first + (fluxSplit[TOP].u - fluxSplit[BOTTOM].u)/cellSizes.second)*dt;
			const double eFlux = ((fluxSplit[RIGHT].e - fluxSplit[LEFT].e)/cellSizes.first + (fluxSplit[TOP].e - fluxSplit[BOTTOM].e)/cellSizes.second)*dt;
			rho.flux(xIdx,yIdx) = rhoFlux;
			v.flux(xIdx,yIdx) = vFlux;
			u.flux(xIdx,yIdx) = uFlux;
			E.flux(xIdx,yIdx) = eFlux;

			// Extra term to compensate for the fact that the cell sizes are not uniform because we are in a cylindrical coordinate system. In Florian (2017), this is the term *** Hr ***.
			EulerContinuity hr;
			const double density = rho.currentTimeStep.GetAt(currentCell);
			const double xVel = u.currentTimeStep.GetAt(currentCell);
			const double yVel = v.currentTimeStep.GetAt(currentCell);
			const double yc = localCellCenterPositions[1].at(currentCell.y);
			const double enthalpy = H.currentTimeStep.GetAt(currentCell);
			
			rho.flux(xIdx,yIdx) += density * yVel / yc * dt;
			u.flux(xIdx,yIdx) = density * yVel * xVel / yc * dt;
			v.flux(xIdx,yIdx) =  density * yVel * yVel / yc * dt;
			E.flux(xIdx,yIdx) =  density * yVel * enthalpy / yc * dt;
		}
	}	
}

void Domain::EmptyFlowDeltaBuffer()
{
	rho.flux.SetAllToValue(0);
	u.flux.SetAllToValue(0);
	v.flux.SetAllToValue(0);
	E.flux.SetAllToValue(0);
	// Note that clearing the others (non-state variables) is not necessary, as they are not written to.
#ifdef _DEBUG
	p.flux.SetAllToValue(0);
	T.flux.SetAllToValue(0);
	H.flux.SetAllToValue(0);
#endif
}

void Domain::SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(const int currentRungeKuttaIter)
{
#ifdef _DEBUG
	std::cout << "Setting next time step values based on runge-kutta and delta buffers for domain '" << name << "'" << std::endl;
	
	// Fail if the buffers are empty; they clearly must be set.
	assert(!rho.rungeKuttaBuffer.IsFilledWithZeroes());
	// assert(!u.rungeKuttaBuffer.IsFilledWithZeroes()); // This can be zero for the initial condition lol.
	// assert(!v.rungeKuttaBuffer.IsFilledWithZeroes()); // This can be zero for the initial condition lol.
	assert(!E.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(rho.flux.IsFilledWithZeroes());
	assert(u.flux.IsFilledWithZeroes());
	assert(v.flux.IsFilledWithZeroes());
	assert(E.flux.IsFilledWithZeroes());
#endif

	const double rkK = 1./(simCase->solverSettings.rungeKuttaOrder-currentRungeKuttaIter); // Runge-kutta factor
	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			// These are state variables, and are explicitly expressed. combine the values in the delta buffers, and apply runge-kutta scaling.
			rho.nextTimeStepBuffer(xIdx, yIdx) = rho.currentTimeStep.GetAt(xIdx, yIdx) - rkK * rho.flux.GetAt(xIdx, yIdx);
			u.nextTimeStepBuffer(xIdx, yIdx) = (rho.currentTimeStep.GetAt(xIdx, yIdx) * u.currentTimeStep.GetAt(xIdx,yIdx) - rkK * u.flux.GetAt(xIdx, yIdx)) / rho.nextTimeStepBuffer(xIdx, yIdx);
			v.nextTimeStepBuffer(xIdx, yIdx) = (rho.currentTimeStep.GetAt(xIdx, yIdx) * v.currentTimeStep.GetAt(xIdx,yIdx) - rkK * v.flux.GetAt(xIdx, yIdx)) / rho.nextTimeStepBuffer(xIdx, yIdx);
			E.nextTimeStepBuffer(xIdx, yIdx) = E.currentTimeStep.GetAt(xIdx, yIdx) - rkK * E.flux.GetAt(xIdx, yIdx);

			// The others are not state variables; they can be calculated using the known variables. Calculate them now.
			//todo: set a build mode where these are not calculated to speed things up, as this is only really necessary for data export.
			p.nextTimeStepBuffer(xIdx, yIdx) = (SpecificHeatRatio()-1) * E.nextTimeStepBuffer.GetAt(xIdx, yIdx) - 0.5*rho.nextTimeStepBuffer(xIdx,yIdx)*pow(u.nextTimeStepBuffer.GetAt(xIdx,yIdx) + v.nextTimeStepBuffer.GetAt(xIdx,yIdx), 2);
			T.nextTimeStepBuffer(xIdx, yIdx) = p.nextTimeStepBuffer.GetAt(xIdx, yIdx) / (GasConstant() * rho.nextTimeStepBuffer.GetAt(xIdx,yIdx));
			H.nextTimeStepBuffer(xIdx,yIdx) = (E.nextTimeStepBuffer.GetAt(xIdx,yIdx) + p.nextTimeStepBuffer(xIdx,yIdx))/rho.nextTimeStepBuffer.GetAt(xIdx,yIdx);
		}
	}


}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y/R!");
	}
}

void Domain::CacheCellSizes()
{
	// Apologies if this function is a little hard to wrap your head around, but this saves a lot of performance pain and repeated code.
	// Since the different axis are independent of eachother, they can be determined individually first. Then, it can be saved for the entire thing.
	cellLengths[0] = std::vector<double>(amountOfCells[0], 0);
	cellLengths[1] = std::vector<double>(amountOfCells[1], 0);
	localCellCenterPositions[0] = std::vector<double>(amountOfCells[0], 0);
	localCellCenterPositions[1] = std::vector<double>(amountOfCells[1], 0);
	// do it for both axes, which are iterates as meshSpacing[axis]
	for (int axis = 0; axis < 2; axis++)
	{
		// As the spacing in the y-direction is not dependent on the x position and vice versa, they can be pre-calculated, and later just populated.
		double previousPosition = 0;
		for (int i = 0; i < amountOfCells[axis]; i++)
		{
			// Get the length of the current cell
			double cellLength = 0;
			switch (meshSpacing[axis].spacingType)
			{
			case EMeshSpacingType::CONSTANT:
				cellLength = meshSpacing[axis].left;
				break;
			case EMeshSpacingType::LINEAR:
				break;

			default:
				throw std::logic_error("Populating domain dimensions is not implemented for this mesh spacing type.");
			}
			
			cellLengths[axis][i] = cellLength;
			// Note that centerPositions[i-1] cannot be used here, because this stores the center positions. On top of it, it would not be defined for the first iteration.
			localCellCenterPositions[axis][i] = previousPosition + 0.5 * cellLength;
			previousPosition += cellLength;
		}
	}
}

void Domain::PopulateSlipConditionGhostCells(const EFace boundary)
{
	const CellIndex ghostOrigin = GetOriginIndexOfBoundary(boundary);
	const auto extent = GetGhostDimensions(boundary);
	const CellIndex domainOrigin = CellIndex(0,0,BOTTOM); // just a small const

	/* Example of what it looks like on left face; Note that the 'Ghost' reference frame points out, relative to the face it is on.
	 *								up	/\ 
	 *		"Domain" reference frame	+ --> X
	 *			
	 *	(-2,Y)	(-1,Y)	|	(0,Y)	(1,Y)	(2,Y)
	 *	(-2,2)	(-1,2)	|	(0,2)	(1,2)	(2,2)	
	 *	(-2,1)	(-1,1)	|	(0,1)	(1,1)	(2,1)
	 *	(-2,0)	(-1,0)	|	(0,0)	(1,0)	(2,0)
	 *
	 *											/\ X
	 *		"Ghost" reference frame		up	<--	+			// TODO: THIS IS MIRRORED!!! FIX.
	 *		
	 *	(2,X)	(1,X)	|	(0,X)	(-1,X)	(-Y,X)
	 *	(2,2)	(1,2)	|	(0,2)	(-1,2)	(-Y,2)
	 *	(2,1)	(1,1)	|	(0,1)	(-1,1)	(-Y,1)
	 *	(2,0)	(1,0)	|	(0,0)	(-1,0)	(-Y,0)
	 */

	// Generally constant in the local x-direction, so y outer loop.
	for (int yLocalIdx = 0; yLocalIdx < extent.second; yLocalIdx++)
	{
		// If it's the top or right face, reverse the direction as that is how the coordinate reference frames are defined. See the drawing below.
		/*
		*				    TOP
		*					/\
		*					|
		*		 + -- > --- > --- > --- +
		*	L	 |			 			|		R
		*	E	/\			 			/\		I
		*	F	 | ->					|  ->	G
		*	T	/\          /\			/\		H
		*		 |          |			|		T
		*		 + -- > --- > --- > --- +
		*				  BOTTOM
		*/
		int yFromGhostCellAwayFromBoundary = yLocalIdx + 1;
		int yOfSourceCellAwayFromBoundary = -yLocalIdx;
		if (boundary == BOTTOM || boundary == LEFT)
		{
			yFromGhostCellAwayFromBoundary *= -1;
			yOfSourceCellAwayFromBoundary *= -1;
		}
		
		for (int xLocalIdx = 0; xLocalIdx < extent.first; xLocalIdx++)
		{
			const CellIndex ghostCellInGhostReferenceFrame = {xLocalIdx, yFromGhostCellAwayFromBoundary, boundary}; // First define it what it is in its local reference frame. Then determine what the corresponding position is relative to the 'origin' of the entire domain. Note that for the ghost cells, the relative location in the reference frame relative to the boundary is negative, and needs to be offset by -1 as well, as 0 is the first positive cell, and not the actual zero-line.
			const CellIndex sourceCellInGhostReferenceFrame =  {xLocalIdx, yOfSourceCellAwayFromBoundary, boundary};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostCellInGhostReferenceFrame, ghostOrigin, domainOrigin );
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem(sourceCellInGhostReferenceFrame, ghostOrigin, domainOrigin);
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(ValidateCellIndex(sourceIndex, false));
#endif
			// Coordinates are expressed in its actual local coordinate frame, so GetAt can be used, and GetReferenceIncludingGhostCells is not necessary.
			rho.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)	=	rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		= - u.rungeKuttaBuffer.GetAt(sourceIndex); // Flipped!
			v.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=   v.rungeKuttaBuffer.GetAt(sourceIndex);
			H.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	T.rungeKuttaBuffer.GetAt(sourceIndex);
		}
	}
}

void Domain::PopulateNoSlipConditionGhostCells(const EFace boundary)
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
			const CellIndex ghostInBoundaryLocalReferenceFrame = {xLocalIdx, ghostY, Opposite(boundary)};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostInBoundaryLocalReferenceFrame, ghostOrigin, {0,0, TOP} );
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem({xLocalIdx, yLocalIdx, Opposite(boundary)}, ghostOrigin, {0,0, TOP});
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(ValidateCellIndex(sourceIndex, false));
#endif

			rho.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)	=	rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		= - u.rungeKuttaBuffer.GetAt(sourceIndex);	// Flipped!
			v.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		= -	v.rungeKuttaBuffer.GetAt(sourceIndex);	// Flipped!
			H.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	T.rungeKuttaBuffer.GetAt(sourceIndex);
		}
	}
}

void Domain::PopulateConnectedGhostCells(const EFace boundary)
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
	Boundary& boundaryOnThisDomain = boundaries.at(boundary);

	// Validate that this is in fact a connected boundary
	if (boundaryOnThisDomain.boundaryType != EBoundaryCondition::CONNECTED || boundaryOnThisDomain.connectedBoundary == nullptr)
		throw std::logic_error("This boundary has not been connected properly.");

#ifdef _DEBUG
	// Also check the opposite boundary
	Boundary* otherBoundary = boundaryOnThisDomain.connectedBoundary;
	if (otherBoundary->boundaryType != EBoundaryCondition::CONNECTED || otherBoundary->connectedBoundary != &boundaryOnThisDomain)
		throw std::logic_error("The boundary is only properly connected one way!");
#endif

	const Domain* otherDomain = boundaryOnThisDomain.connectedBoundary->domain;
	const bool bVertical = (boundary == LEFT || boundary == RIGHT);

#ifdef _DEBUG
	// Check if they have the same amount of cells.
	assert(amountOfCells[bVertical] == otherDomain->amountOfCells[bVertical]);
#endif

	CellIndex ghostOrigin = GetOriginIndexOfBoundary(boundary);
	CellIndex complementOrigin = otherDomain->GetOriginIndexOfBoundary(Opposite(boundary));
	auto extent = GetGhostDimensions(boundary);
	const CellIndex DomainOrigin = CellIndex(0,0,BOTTOM); // just a small const

	// Generally constant in the local x-direction, so y outer loop.
	for (int yLocalIdx = 0; yLocalIdx < extent.second; yLocalIdx++)
	{
		// If it's the top or right face, reverse the direction as that is how the coordinate reference frames are defined. See the drawing below.
		/*
		*				    TOP
		*					/\
		*					|
		*		 + -- > --- > --- > --- +
		*	L	 |			 			|		R
		*	E	/\			 			/\		I
		*	F	 | ->					|  ->	G
		*	T	/\          /\			/\		H
		*		 |          |			|		T
		*		 + -- > --- > --- > --- +
		*				  BOTTOM
		*/
		
		int yFromGhostCellAwayFromBoundary = yLocalIdx + 1;
		int yOfSourceCellAwayFromBoundary = -yLocalIdx;
		if (boundary == BOTTOM || boundary == LEFT)
		{
			yFromGhostCellAwayFromBoundary *= -1;
			yOfSourceCellAwayFromBoundary *= -1;
		}
		
		for (int xLocalIdx = 0; xLocalIdx < extent.first; xLocalIdx++)
		{
			const CellIndex ghostCellInGhostReferenceFrame = {xLocalIdx, yFromGhostCellAwayFromBoundary, boundary}; // First define it what it is in its local reference frame. Then determine what the corresponding position is relative to the 'origin' of the entire domain. Note that for the ghost cells, the relative location in the reference frame relative to the boundary is negative, and needs to be offset by -1 as well, as 0 is the first positive cell, and not the actual zero-line.
			const CellIndex sourceCellInOtherBoundaryReferenceFrame =  {xLocalIdx, -yOfSourceCellAwayFromBoundary, Opposite(boundary)};
			const CellIndex ghostIndex = TransformToOtherCoordinateSystem(ghostCellInGhostReferenceFrame, ghostOrigin, DomainOrigin );
			const CellIndex sourceIndex = TransformToOtherCoordinateSystem(sourceCellInOtherBoundaryReferenceFrame, complementOrigin, DomainOrigin); // Defined relative to the other domain's (0,0)
			
#ifdef _DEBUG
			assert(ValidateCellIndex(ghostIndex, true));
			assert(otherDomain->ValidateCellIndex(sourceIndex, false));
#endif

			rho.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	otherDomain->rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	otherDomain->p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=   otherDomain->u.rungeKuttaBuffer.GetAt(sourceIndex);	
			v.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=  	otherDomain->v.rungeKuttaBuffer.GetAt(sourceIndex);	
			H.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	otherDomain->H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	otherDomain->E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.GetReferenceIncludingGhostCells(ghostIndex)		=	otherDomain->T.rungeKuttaBuffer.GetAt(sourceIndex);
			
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

bool Domain::AreAllBoundariesSet() const
{
	if (boundaries.size() < 4 || boundaries.size() > 4)
		return false;
	for (const auto& it : boundaries)
	{
		const Boundary& boundary = it.second;
		// This could be moved into boundary.h/boundary.cpp
		if (boundary.boundaryType == EBoundaryCondition::NOT_SET)
			return false;
		if (boundary.boundaryType == EBoundaryCondition::CONNECTED)
		{
			// Ensure the other one is also set as connected
			if (boundary.connectedBoundary->boundaryType != EBoundaryCondition::CONNECTED)
				return false;
		}
	}
	return true;
}

std::pair<int, int> Domain::GetGhostDimensions(EFace boundary)
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

