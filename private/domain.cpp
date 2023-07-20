#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cassert>
#include <ios>
#include <iostream>
#include <cmath>

#include "cell_values_container.h"
#include "euler_container.h"
#include "flux_splitting.h"


Domain::Domain(const std::string& name, SimCase* simCase, const Position& position, const std::pair<double, double> sizeArg, const std::pair<MeshSpacingSolution, MeshSpacingSolution> meshSpacingArg, const EInitialisationMethod initialisationMethod, const int ghostCellDepth) :
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
	meshSpacing[0] = MeshSpacingSolution(meshSpacingArg.first);
	meshSpacing[1] = MeshSpacingSolution(meshSpacingArg.second);

	for (auto& conservationEq : eulerConservationTerms)
		conservationEq.Resize(amountOfCells[0], amountOfCells[1]);

	CacheCellSizes();

}

CellIndex Domain::InvertPositionToIndex(const Position pos, Position& distanceFromCenterOut) const
{
#ifdef _DEBUG
	if (pos.x < 0 || pos.x > size[0] || pos.y < 0 || pos.y > size[1])
		throw std::logic_error("Cannot invert position that is outside of the bounds of the domain.");
#endif
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
		if (IsCloserToLeftThanToRight(pos.y, localCellCenterPositions[1].at(leftFromIndex), localCellCenterPositions[1].at(leftFromIndex+1 )))
		{
			out.y = out.x = static_cast<int>(leftFromIndex);
			distanceFromCenterOut.y = pos.y - localCellCenterPositions[1].at(leftFromIndex);
		}
		else
		{
			out.y = static_cast<int>(leftFromIndex) + 1;
			distanceFromCenterOut.y = pos.y - localCellCenterPositions[1].at(leftFromIndex + 1);
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
	double volume = cellDX * M_PI * (std::pow(rOuter, 2) - std::pow(rInner, 2));
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
	const double ESet = pSet / (gamma - 1.0) + 0.5 * rhoSet * (std::pow(uSet, 2) + std::pow(vSet, 2));
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
		const EFace face = it.first;
		const Boundary& boundary = it.second;

#ifdef _DEBUG
		std::cout << "Updating Ghost Cells for domain '" << name << "', side" << FaceToString(face) << std::endl;
#endif

		// Note that the below functions all work in relative coordinate frames.
		// Technically a performance gain could be achieved by not calling a transformation on the frame every time, instead directly accessing those entries directly by looping. However, this is more legible.
		switch (boundary.boundaryType) {
		case EBoundaryCondition::NOT_SET:
				throw std::logic_error("Cannot update Ghost cell if boundary condition is not set.");
			case EBoundaryCondition::SLIP:
				PopulateSlipConditionGhostCells(face);
				break;
			case EBoundaryCondition::NO_SLIP:
				PopulateNoSlipConditionGhostCells(face);
				break;
			case EBoundaryCondition::CONNECTED:
				PopulateConnectedGhostCells(face);
				break;
		default:
				throw std::logic_error("Updating ghost cells is not implemented for this type of boundary.");
		}
	}
}

void Domain::CacheEulerConservationTerms(const double dt)
{
	// Depending on whether or not there are shock fronts, the integration scheme changes. Determining whether or not there are shock fronts is done by looking at the MUSCL interpolated values of the field quantities. So, first calculate & cache the MUSCL vales.
	const SolverSettings& solverSettings = simCase->solverSettings;

#ifdef _DEBUG
	std::cout << "Calculating flow delta, and populating Flow Delta buffer of domain '" << name << "'" << std::endl;
	
	// Ensure that the buffers are actually empty to start with
	if (!eulerConservationTerms[0].IsFilledWithZeroes())
		throw std::logic_error("conservation of mass buffer is non-empty");
	if (!eulerConservationTerms[1].IsFilledWithZeroes())
		throw std::logic_error("conservation of x-momentum buffer is non-empty");
	if (!eulerConservationTerms[2].IsFilledWithZeroes())
		throw std::logic_error("conservation of y-momentum buffer is non-empty");
	if (!eulerConservationTerms[2].IsFilledWithZeroes())
		throw std::logic_error("conservation of energy buffer is non-empty");
#endif

	// Do flux splitting on all the faces. Use hanel if flow is sonic in either cell, ausm_dv if not.
	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			const CellIndex cix(xIdx, yIdx);
			double gamma = SpecificHeatRatio();
			MUSCLBuffer rhoMUSCL, uMUSCL, vMUSCL, pMUSCL, eMUSCL, hMUSCL;
			constexpr EFace faces[4] = {LEFT, RIGHT, TOP, BOTTOM}; // Decaying the faces to integers makes the calcuations more human readable. indexes for sides ---- 0: right; 1: left; 2: top; 3: down; 
			constexpr EAxisDirection faceNormalDirections[2] = {EAxisDirection::POSITIVE, EAxisDirection::NEGATIVE};

			// Getting the MUSCL interpolations for the values at this cell's faces.
			// This could be made quicker with caching, but at this makes the program significantly more complicated. It has therefore been decided to keep it like this.
			for (const auto f : faces)
			{
				for (const auto n : faceNormalDirections)
				{
					rhoMUSCL.At(f, n) = rho.rungeKuttaBuffer.GetMUSCLInterpolationForFace(cix, f, n, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
					uMUSCL.At(f, n) = u.rungeKuttaBuffer.GetMUSCLInterpolationForFace(cix, f, n, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
					vMUSCL.At(f, n) = v.rungeKuttaBuffer.GetMUSCLInterpolationForFace(cix, f, n, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
					pMUSCL.At(f, n) = p.rungeKuttaBuffer.GetMUSCLInterpolationForFace(cix, f, n, solverSettings.MUSCLBias, solverSettings.fluxLimiterType);
					// rho, u, v and p's MUSCL values determined based on their state variable (field). the other variables are not, but they can be calculated directly from the state values in those buffers, as they are non-state variables. So, do this for energy and enthalpy.
					eMUSCL.At(f, n) = pMUSCL.GetAt(f, n) / (gamma -1) +  0.5*rhoMUSCL.GetAt(f, n) * (std::pow(uMUSCL.GetAt(f, n), 2) + std::pow(vMUSCL.GetAt(f, n), 2));
					hMUSCL.At(f, n) = (eMUSCL.GetAt(f, n) + pMUSCL.GetAt(f, n))  / rhoMUSCL.GetAt(f, n); 
				}
			}
			
			// todo: Florian (2017) did this by setting a global flag. If it had a shockwave in one direction, all of the faces were marked, as well as the cell it neighbours. The approach currently employed only sets this face, and not the neighbours. Is this admissable?
			bool shockwavePresent = false;
			for (const auto f : faces)
			{
				const double cNeg = sqrt(SpecificHeatRatio() * pMUSCL.At(f, EAxisDirection::NEGATIVE));
				const double cPos = sqrt(SpecificHeatRatio() * pMUSCL.At(f, EAxisDirection::POSITIVE));
				// using !=, Logical xor operator: one of them is supersonic, the other is not. More intuitive than ( (uL-cL)>0 && (uR-cR)<0 ) || ( (uL+cL)>0 && (uR+cR)<0 )
				if (f == LEFT || f == RIGHT)
					if ((std::abs(uMUSCL.GetAt(f, EAxisDirection::POSITIVE)) > cPos) != (std::abs(uMUSCL.GetAt(f, EAxisDirection::NEGATIVE)) > cNeg))
						shockwavePresent = true;
				if (f == TOP || f == BOTTOM)
					if ((std::abs(vMUSCL.GetAt(f, EAxisDirection::POSITIVE)) > cPos) != (std::abs(vMUSCL.GetAt(f, EAxisDirection::NEGATIVE)) > cNeg))
						shockwavePresent = true;
			}
			
			// Now, for every face, determine which flux scheme to use based on whether or not it is critical (sonic), and perform the flux transfer.
			EulerContinuity continuityAtFace[4]; // These here contain 4 terms for the euler equations. They all represent the flux of those variables at each side of the cell.
			for (const auto f : faces)
			{
				// Get the value of the flow going in the positive normal direction and negative normal direction at this specific face. 
				CellValues positiveNormalFlow;
				positiveNormalFlow.density = rhoMUSCL.At(f, EAxisDirection::POSITIVE);
				positiveNormalFlow.u = uMUSCL.At(f, EAxisDirection::POSITIVE);
				positiveNormalFlow.v = vMUSCL.At(f, EAxisDirection::POSITIVE);
				positiveNormalFlow.e = eMUSCL.At(f, EAxisDirection::POSITIVE);
				positiveNormalFlow.h = hMUSCL.At(f, EAxisDirection::POSITIVE);
				positiveNormalFlow.p = pMUSCL.At(f, EAxisDirection::POSITIVE);
				
				CellValues negativeNormalFlow;
				negativeNormalFlow.density = rhoMUSCL.At(f, EAxisDirection::NEGATIVE);
				negativeNormalFlow.u = uMUSCL.At(f, EAxisDirection::NEGATIVE);
				negativeNormalFlow.v = vMUSCL.At(f, EAxisDirection::NEGATIVE);
				negativeNormalFlow.e = eMUSCL.At(f, EAxisDirection::NEGATIVE);
				negativeNormalFlow.h = hMUSCL.At(f, EAxisDirection::NEGATIVE);
				negativeNormalFlow.p = pMUSCL.At(f, EAxisDirection::NEGATIVE);

				// Set the flux split for this side (declared above in the array).
				// If there's a shockwave, use Hanel. If there is no shockwave, use AUSM DV.
				const bool bIsVertical = (f == TOP || f == BOTTOM);
				if (shockwavePresent)
					continuityAtFace[f] = HanelFluxSplitting(negativeNormalFlow, positiveNormalFlow, bIsVertical, gamma, solverSettings.entropyFix);
				else // No shockwave
					continuityAtFace[f] = AUSMDVFluxSplitting(negativeNormalFlow, positiveNormalFlow, bIsVertical,gamma, solverSettings.AUSMSwitchBias, solverSettings.entropyFix);

			}

			// Total accumulation is what goes in - what goes out
			const auto cellSizes = GetCellSizes(cix);
			EulerContinuity accumulationAuto = ((continuityAtFace[RIGHT] - continuityAtFace[LEFT])/cellSizes.first + (continuityAtFace[TOP] - continuityAtFace[BOTTOM])/cellSizes.second)*dt; 
			EulerContinuity accumulation;
			accumulation.mass = ((continuityAtFace[RIGHT].mass - continuityAtFace[LEFT].mass)/cellSizes.first + (continuityAtFace[TOP].mass - continuityAtFace[BOTTOM].mass)/cellSizes.second)*dt;
			accumulation.momentumX = ((continuityAtFace[RIGHT].momentumX - continuityAtFace[LEFT].momentumX)/cellSizes.first + (continuityAtFace[TOP].momentumX - continuityAtFace[BOTTOM].momentumX)/cellSizes.second)*dt;
			accumulation.momentumY = ((continuityAtFace[RIGHT].momentumY - continuityAtFace[LEFT].momentumY)/cellSizes.first + (continuityAtFace[TOP].momentumY - continuityAtFace[BOTTOM].momentumY)/cellSizes.second)*dt;
			accumulation.energy = ((continuityAtFace[RIGHT].energy - continuityAtFace[LEFT].energy)/cellSizes.first + (continuityAtFace[TOP].energy - continuityAtFace[BOTTOM].energy)/cellSizes.second)*dt;
			// dy = dr in this case, if you compensate for the squashification.

			// Extra term to compensate for the fact that the cell sizes are not uniform because we are in a cylindrical coordinate system. In Florian (2017), this is the term *** Hr ***.
			const double density = rho.currentTimeStep.GetAt(cix);
			const double xVel = u.currentTimeStep.GetAt(cix);
			const double yVel = v.currentTimeStep.GetAt(cix);
			const double yc = localCellCenterPositions[1].at(cix.y);
			const double enthalpy = H.currentTimeStep.GetAt(cix);
			accumulation.mass += density * yVel / yc * dt;
			accumulation.momentumX += density * yVel * xVel / yc * dt;
			accumulation.momentumY +=  density * yVel * yVel / yc * dt;
			accumulation.energy +=  density * yVel * enthalpy / yc * dt;

			eulerConservationTerms[0](xIdx, yIdx) = accumulation.mass;
			eulerConservationTerms[1](xIdx, yIdx) = accumulation.momentumX;
			eulerConservationTerms[2](xIdx, yIdx) = accumulation.momentumY;
			eulerConservationTerms[3](xIdx, yIdx) = accumulation.energy;
		}
	}	
}

void Domain::EmptyEulerConservationCache()
{
	for (auto& buf : eulerConservationTerms)
		buf.SetAllToValue(0);
}

void Domain::SetNextTimeStepValuesBasedOnCachedEulerContinuities(const int currentRungeKuttaIter)
{
#ifdef _DEBUG
	std::cout << "Setting next time step values based on runge-kutta and delta buffers for domain '" << name << "'" << std::endl;
	
	// Fail if the buffers are empty; they clearly must be set.
	assert(!rho.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(!p.rungeKuttaBuffer.IsFilledWithZeroes());
	assert(!E.rungeKuttaBuffer.IsFilledWithZeroes());
#endif

	const double rkK = 1./(simCase->solverSettings.rungeKuttaOrder-currentRungeKuttaIter); // Runge-kutta factor
	for (int xIdx = 0; xIdx < amountOfCells[0]; xIdx++)
	{
		for (int yIdx = 0; yIdx < amountOfCells[1]; yIdx++)
		{
			const CellIndex cix = {xIdx, yIdx};
			const EulerContinuity accumulation = {eulerConservationTerms[0].GetAt(cix), eulerConservationTerms[1].GetAt(cix), eulerConservationTerms[2].GetAt(cix), eulerConservationTerms[3].GetAt(cix)};
			// Convert the conservation equations back into actual variables here.
			// These are state variables, and are explicitly expressed. combine the values in the delta buffers, and apply runge-kutta scaling.
			const double deltaRho = -rkK * accumulation.mass;
			const double deltaE = -rkK * accumulation.momentumY;
			
			const double nextRho = rho.currentTimeStep.GetAt(cix) + deltaRho;
			const double nextU = (rho.currentTimeStep.GetAt(cix) * u.currentTimeStep.GetAt(cix) - rkK * accumulation.momentumX) / nextRho;
			const double nextV = (rho.currentTimeStep.GetAt(cix) * v.currentTimeStep.GetAt(cix) - rkK * accumulation.momentumY) / nextRho;
			const double nextE = E.currentTimeStep.GetAt(cix) + deltaE;

			rho.nextTimeStepBuffer(cix) = nextRho;
			u.nextTimeStepBuffer(cix) = nextU;
			v.nextTimeStepBuffer(cix) = nextV;
			E.nextTimeStepBuffer(cix) = nextE;
			
			// The others are not state variables; they can be calculated using the known variables. Calculate them now.
			//todo: set a build mode where these are not calculated unless a record has been set.
			const double nextP = (SpecificHeatRatio()-1) * (nextE - 0.5*nextRho*(nextU*nextU + nextV*nextV));
			p.nextTimeStepBuffer(cix) = nextP;
			const double nextT = nextP / (GasConstant() * nextRho);
			T.nextTimeStepBuffer(cix) = nextT;
			const double nextH = (nextE + nextP)/nextRho;
			H.nextTimeStepBuffer(cix) = nextH;

#ifdef _DEBUG
			assert(rho.nextTimeStepBuffer(cix) > 0);
			assert(E.nextTimeStepBuffer(cix) > 0);
			assert(p.nextTimeStepBuffer(cix) > 0);
			assert(H.nextTimeStepBuffer(cix) > 0);
			assert(T.nextTimeStepBuffer(cix) > 0); // Kelvins, not C.
#endif
			
		}
	}
}

void ValidateAxisInput(const int axis)
{
	if (axis > 1 || axis < 0)
	{
		throw std::invalid_argument("Cannot get access axis with the required index; there are only 2, X and Y (R)!");
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
			rho.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)	=	rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		= - u.rungeKuttaBuffer.GetAt(sourceIndex); // Flipped!
			v.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=   v.rungeKuttaBuffer.GetAt(sourceIndex);
			H.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	T.rungeKuttaBuffer.GetAt(sourceIndex);
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

			rho.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)	=	rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		= - u.rungeKuttaBuffer.GetAt(sourceIndex);	// Flipped!
			v.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		= -	v.rungeKuttaBuffer.GetAt(sourceIndex);	// Flipped!
			H.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	T.rungeKuttaBuffer.GetAt(sourceIndex);
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

			rho.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	otherDomain->rho.rungeKuttaBuffer.GetAt(sourceIndex);
			p.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	otherDomain->p.rungeKuttaBuffer.GetAt(sourceIndex);
			u.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=   otherDomain->u.rungeKuttaBuffer.GetAt(sourceIndex);	
			v.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=  	otherDomain->v.rungeKuttaBuffer.GetAt(sourceIndex);	
			H.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	otherDomain->H.rungeKuttaBuffer.GetAt(sourceIndex);
			E.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	otherDomain->E.rungeKuttaBuffer.GetAt(sourceIndex);
			T.rungeKuttaBuffer.AtWithGhostCells(ghostIndex)		=	otherDomain->T.rungeKuttaBuffer.GetAt(sourceIndex);
			
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
	std::map<EFace, bool> truth;
	
	if (boundaries.size() < 4)
		return false;
	if (boundaries.size() > 4)
		throw std::logic_error("Too many boundaries defined! This should not be possible.");
	
	for (const auto& it : boundaries)
	{
		if (truth.count(it.first))
			throw std::logic_error("The side" + FaceToString(it.first) + "Already exists!");
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
	case TOP:
		return {amountOfCells[0],nGhost};
	case RIGHT:
		return {amountOfCells[1],nGhost};
	default:
		throw std::logic_error("GetGhostOrigin is not implemented for this boundary location.");
	}
}

