#include "reed_valve.h"

#include <cassert>
#include <numeric>

#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cmath>
#include "fem_deformation.h"

#include "AuxFunctions.h"

#define HOLE_FACTOR 0.9

#define DAMPING_C1 5.0E-8 //5.0E-7  
#define DAMPING_C2 0.0E-8 //2.0E-8 // Why on earth is this 0??
#define DAMPING_C3 0.0007
#define SAMPLING_DEPTH_FOR_OUT_OF_DOMAIN 2

ReedValve::ReedValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile, const bool bMirrored) :
	IValve(intoDomain, outOfDomain, boundary, positionAlongBoundary),
	bMirrored(bMirrored),
	amountOfFixedNodes(amountOfFixedNodes),
	amountOfFreeNodes(amountOfFreeSections),
	lengthOfFreeSection(lengthOfFreeSection),
	lengthOfFixedSections(lengthOfFixedSections),
	beamProfile_(beamProfile)
{
	positionMirrorModifier_ = bMirrored ? -1 : 1;
	hingePositionInDomain = intoDomain->PositionAlongBoundaryToCoordinate(boundary, positionAlongBoundary, 0);
	holeEndPositionAlongBoundary = positionAlongBoundary + (lengthOfFixedSections + lengthOfFreeSection) * positionMirrorModifier_;
	holeEndPositionInDomain = intoDomain->PositionAlongBoundaryToCoordinate(boundary, holeEndPositionAlongBoundary , 0);

	hingePositionIndex_ = intoDomain->InvertPositionToIndex(hingePositionInDomain);
	holeEndPositionIndex_ = intoDomain->InvertPositionToIndex(holeEndPositionInDomain);

	// Instead, only create it at Register().
	//fem_ = FemDeformation(amountOfFreeSections, amountOfFixedNodes, beamProfile, lengthOfFreeSection, lengthOfFixedSections, intoDomain->simCase->dt, Opposite(boundary));
}

void ReedValve::CalculateForceOnNodesFromPressure(std::vector<double>& forceVectorOut, const EFieldQuantityBuffer bufferName) const
{
	// Only do this for the nodes that are considered 'free'.
	// Note that the beam connecting the last fixed and the first free node is still considered 'fixed', it cannot create loading, as this would be impossible to distribute between the two nodes.
	/*
	 * Graphical representation:
	 *   S1  S2  S3  S4  S5
	 *               \/  \/
	 * X---X---X---O---O---O
	 * N1  N2  N3  N4  N5  N6
	 *
	 * f on N3 = 0.5 * S3
	 */
#ifdef _DEBUG
	// Reading from a pressure field, so pressure cannot be 0 (possibly unset)
	if (intoDomain_->p.bufferMap.at(bufferName).IsFilledWithZeroes())
		throw std::logic_error("intoDomain's pressure field is not initialised.");
	if (outOfDomain_->p.bufferMap.at(bufferName).IsFilledWithZeroes())
		throw std::logic_error("outOfDomain's pressure field is not initialised.");
#endif
	
	forceVectorOut.resize(fem_.amountOfNodes * N_DOF_PER_NODE);
	
	// This iterates over the beam section elements, but we need the indices to determine the positions. Hence, up to amountOfNodes-1
	for (int nodeIdx = fem_.fixedNodes; nodeIdx < fem_.amountOfNodes - 1; nodeIdx++)
	{
		// Sample the pressures where the element is in the physical domain. To get the position where the beam section is in the total domain, the positions of the nodes that it spans between are averaged, it is converted to the reference frame of the domain (instead of the reed valve), and then added to the actual position of the valve in the domain.
		const Position beamSectionCenterPositionLocal = (fem_.nodePositionsRelativeToRoot[nodeIdx] + fem_.nodePositionsRelativeToRoot[nodeIdx + 1]) * 0.5;
		const Position beamSectionCenterPositionInDomain = TransformToOtherCoordinateSystem(beamSectionCenterPositionLocal, hingePositionInDomain, {0,0});
		const double pressureAtBeamCenter = intoDomain_->p.GetInterpolatedValueAtPosition(beamSectionCenterPositionInDomain, bufferName);

		// The pressure is now sampled in a very similar manner in the the domain that this valve source from.
		// To get the position, the fact that the position of the node is known in a local coordinate system is used to essentially 'mirror' it over the boundary, into the other domain.
		const double depthIntoOtherDomain = 2; // Right now this is just a fixed number. // todo make this a parameter.
		Position ambientSamplePositionLocal;
		switch (boundary_)
		{
		case LEFT:
			ambientSamplePositionLocal = {outOfDomain_->size[0] - depthIntoOtherDomain, beamSectionCenterPositionInDomain.y};
			break;
		case RIGHT:
			ambientSamplePositionLocal = {depthIntoOtherDomain, beamSectionCenterPositionInDomain.y};
			break;
		case TOP:
			ambientSamplePositionLocal = {beamSectionCenterPositionInDomain.x, depthIntoOtherDomain};
			break;
		case BOTTOM:
			ambientSamplePositionLocal = {beamSectionCenterPositionInDomain.x, outOfDomain_->size[1] - depthIntoOtherDomain};
			break;
		default:
			throw std::logic_error("Sampling position in other domain is not implemented for this boundary type.");
		}

		double pressureAtSink =  outOfDomain_->p.GetInterpolatedValueAtPosition(ambientSamplePositionLocal, bufferName);
		
		const double deltaPressureWithAmbient = pressureAtBeamCenter - pressureAtSink;

		// We're only interested in the (locally) vertical component. hence, determine theta, the angle it makes relative to the (local) horizontal axis.
		const Position deltaPosition = fem_.nodePositionsRelativeToRoot[nodeIdx + 1] - fem_.nodePositionsRelativeToRoot[nodeIdx];
		double cosTheta = deltaPosition.x / deltaPosition.Distance(); // Just really simple pythagoras.
		
		double forceOnElement = deltaPressureWithAmbient * fem_.beamSections.at(nodeIdx).topOrBottomSurfaceArea * cosTheta;

		// The forces is assumed to be equally distributed over the two different nodes.
		// Note that only the local y-force is calculated, so the others are left empty.
		forceVectorOut[nodeIdx * N_DOF_PER_NODE] += 0.5* forceOnElement;
		forceVectorOut[(nodeIdx + 1) * N_DOF_PER_NODE] += 0.5* forceOnElement;
		
	}
}

void ReedValve::OnRegister()
{
	fem_ = FemDeformation(amountOfFreeNodes, amountOfFixedNodes, beamProfile_, lengthOfFreeSection, lengthOfFixedSections, intoDomain_->simCase->dt, boundary_);
	SetSourceCellIndices(sourceCellIndices, boundary_, positionAlongBoundary_, lengthOfFreeSection, lengthOfFixedSections);
	
}

void ReedValve::SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EFace boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections) const
{
	// calculate the 'starting position' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection * (1 - HOLE_FACTOR);
	auto posStart = intoDomain_->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart, 0);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain_->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd, 0);

	CellIndex sourceStartIndexOnBoundary = intoDomain_->InvertPositionToIndex(posStart);		// The index (location) on the boundary where the valve starts creating a source term.
	CellIndex sourceEndIndexOnBoundary = intoDomain_->InvertPositionToIndex(posEnd);		// The index (location) on the boundary where the valve stops creating a source term.

	// Determine all the positions between the two points, and save them as source terms by interpolating a line between the two.
	// This implementation assumes it is either perfectly horizontal or vertical, and does not allow for slanted lines or other profiles.
	bool bHorizontalDifference = (sourceStartIndexOnBoundary.x != sourceEndIndexOnBoundary.x);
	bool bVerticalDifference = (sourceStartIndexOnBoundary.y != sourceEndIndexOnBoundary.y);

	if (bHorizontalDifference == bVerticalDifference)
	{
		// Note that this will also throw if the total size is only 1x1!m Maybe need to make an exception for htis, but generally you wouldn't want this anyway.
		throw std::logic_error("A reed valve cannot have a source term in both directions!");
	}

	// Populate the source cell indices array with the value that does not change as a constant line.
	if (bHorizontalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.x; i < sourceEndIndexOnBoundary.x; i++)
		{
			sourceCellIndicesOut.emplace_back(i, sourceStartIndexOnBoundary.y);
		}
	}
	else if (bVerticalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.y; i < sourceEndIndexOnBoundary.y; i++)
		{
			sourceCellIndicesOut.emplace_back(sourceStartIndexOnBoundary.x, i);
		}
	}

	#ifdef _DEBUG
	if (bHorizontalDifference)
		assert(sourceCellIndices.size() < intoDomain_->size[0]);
	else if (bVerticalDifference)
		assert(sourceCellIndices.size() < intoDomain_->size[1]);
	#endif
}

void ReedValve::Update()
{
	#ifdef _DEBUG
	if (IsCloseToZero(fem_.dt - 0))
		throw std::logic_error("FEM module is not initialised.");
	#endif
	
	// TODO: Right now creates & destroys them every time. Possible optimisation would be to cache them?
	std::vector<double> forcesOnNodes; // The forces on each node in the local vertical (y positive!) direction.
	std::vector<double> newDeflection; // The solution to the FEM will be here. it goes [x_node1, y_node1, x_node2, y_node2, (...)]
	CalculateForceOnNodesFromPressure(forcesOnNodes, RUNGE_KUTTA);
	CalculateAerodynamicDamping(forcesOnNodes);
	fem_.CalculateNewDeflections(newDeflection, forcesOnNodes);
	fem_.UpdatePositions(newDeflection);
}

double ReedValve::GetMassFlowRate() const
{
	double tipDeflection = GetTipDeflection();
	if (IsCloseToZero(tipDeflection) || tipDeflection < 0)
	{
		// This setting is a bit artifical, as the tip can also be arching, still leaving some mass flow. It comes close enough to reality though.
		return 0;
	}
	
	// todo: replace by SetPressureRatios
	double averagePressureIntoDomain = GetAverageFieldQuantityAroundValve(intoDomain_->p, CURRENT_TIME_STEP, true);
	double averagePressureOutOfDomain = GetAverageFieldQuantityAroundValve(outOfDomain_->p, CURRENT_TIME_STEP, true);
	double averageDensityOutOfDomain = GetAverageFieldQuantityAroundValve(outOfDomain_->rho, CURRENT_TIME_STEP);

#ifdef _DEBUG
	if (averagePressureIntoDomain < 0 || averagePressureOutOfDomain < 0 || averageDensityOutOfDomain < 0 )
		throw std::logic_error("Average pressures or densities cannot be lower than 0.");
#endif

	const double pressureRatio = fmax(0.0,fmin(averagePressureIntoDomain/averagePressureOutOfDomain,1.0));
	if (pressureRatio > 1.0)
	{
		// This only makes sense if the reed valve is completely closed. There is a very small moment where this pressure gradient may be this way while it is still open, but this is disregarded.
		return 0;
	}

	const double referenceArea = FukanariReferenceArea();
	const double gamma = intoDomain_->SpecificHeatRatio();
	// If the critical pressure ratio is reached, the flow is considered choked, so the mass flow can no longer be increased!
	const double criticalPressureRatio = pow(2.0/(gamma+1.0),gamma/(gamma-1.0));
	const double dischargeCoefficient = DischargeCoefficient();
	
	if (pressureRatio > criticalPressureRatio) // && pratio <= 1.0), flow is choked.
	{
		double x = dischargeCoefficient*referenceArea*averageDensityOutOfDomain*pow(pressureRatio,1.0/gamma)*sqrt(2.0*gamma*averagePressureOutOfDomain/averageDensityOutOfDomain/(gamma-1.0)*(1.0-pow(pressureRatio,(gamma-1.0)/gamma)));

#ifdef _DEBUG
		assert(!isnan(x));
#endif
		return x;
	}
	else // pressureRatio > 0.0 && pressureRatio <= criticalPressureRatio, flow is not choked.
	{
		double x = dischargeCoefficient*referenceArea*gamma*averagePressureOutOfDomain/sqrt(gamma*averagePressureOutOfDomain/averageDensityOutOfDomain)*pow(2.0/(gamma+1.0),0.5*(gamma+1.0)/(gamma-1.0));

#ifdef _DEBUG
		assert(!isnan(x));
#endif
		return x;
	}
}

void ReedValve::FillBuffer()
{
	// First part: calculating the total mass flow rate. in legacy code, equal to calculate_mfr function.
	double tipDeflection = GetTipDeflection();

	// Early exit with zeroes if the valve is closed.
	if (IsCloseToZero(tipDeflection) || tipDeflection < 0)
	{
		// This setting is a bit artifical, as the tip can also be arching, still leaving some mass flow. It comes close enough to reality though.
		for (auto& a : sourceTermBuffer_)
			a = EulerContinuity();
		return;
	}

	double averagePressureIntoDomain = GetAverageFieldQuantityAroundValve(intoDomain_->p, CURRENT_TIME_STEP, true);
	double averagePressureOutOfDomain = GetAverageFieldQuantityAroundValve(outOfDomain_->p, CURRENT_TIME_STEP, true);
	double averageDensityOutOfDomain = GetAverageFieldQuantityAroundValve(outOfDomain_->rho, CURRENT_TIME_STEP, false);

#ifdef _DEBUG
	if (averagePressureIntoDomain < 0 || averagePressureOutOfDomain < 0 || averageDensityOutOfDomain < 0 )
		throw std::logic_error("Average pressures or densities cannot be lower than 0.");
#endif

	const double averagePressureRatio = fmax(0.0,fmin(averagePressureIntoDomain/averagePressureOutOfDomain,1.0));

	// Early exit that does not allow flow reversal.
	if (averagePressureRatio > 1.0)
	{
		// This only makes sense if the reed valve is completely closed. There is a very small moment where this pressure gradient may be this way while it is still open, but this is disregarded.
		for (auto& a : sourceTermBuffer_)
			a = EulerContinuity();
		return;
	}

	const double referenceArea = FukanariReferenceArea();
	const double gamma = intoDomain_->SpecificHeatRatio();
	// If the critical pressure ratio is reached, the flow is considered choked, so the mass flow can no longer be increased!
	const double criticalPressureRatio = pow(2.0/(gamma+1.0),gamma/(gamma-1.0));
	const double dischargeCoefficient = DischargeCoefficient();

#ifdef _DEBUG
	assert(!IsCloseToZero(averagePressureRatio));
	assert(!IsCloseToZero(criticalPressureRatio));
#endif

	double totalMassFlowRate;
	
	if (averagePressureRatio > criticalPressureRatio) // && pratio <= 1.0), flow is choked.
	{
		totalMassFlowRate = dischargeCoefficient*referenceArea*averageDensityOutOfDomain*pow(averagePressureRatio,1.0/gamma)*sqrt(2.0*gamma*averagePressureOutOfDomain/averageDensityOutOfDomain/(gamma-1.0)*(1.0-pow(averagePressureRatio,(gamma-1.0)/gamma)));

#ifdef _DEBUG
		assert(!isnan(totalMassFlowRate));
#endif
	}
	else // pressureRatio > 0.0 && pressureRatio <= criticalPressureRatio, flow is not choked.
	{
		totalMassFlowRate = dischargeCoefficient*referenceArea*gamma*averagePressureOutOfDomain/sqrt(gamma*averagePressureOutOfDomain/averageDensityOutOfDomain)*pow(2.0/(gamma+1.0),0.5*(gamma+1.0)/(gamma-1.0));

#ifdef _DEBUG
		assert(!isnan(totalMassFlowRate));
#endif
	}

	/*** Here is what old cold has in compute_cell_source ***/
	
#ifdef _DEBUG
	assert(sourceCellsIndices_.size() == sourceTermBuffer_.size());
	assert(intoDomain_ != nullptr);
	if (sourceCellsIndices_.empty())
		throw std::logic_error("No source cells are set for this valve.");

#endif

	auto volumes = std::vector<double>(sourceTermBuffer_.size(), 0);
	//todo: cache this, if the sources cannot be moved.

	// Assume that the mass flow is equally distributed over all the cells, weighted by cell volume.
	// todo: weighing can be done relative to the tip deflection. Would have to be validated though.
	// todo: separate function for calculating cell volumes.
	for (size_t i = 0; i < sourceTermBuffer_.size(); i++)
	{
		const CellIndex& cix = sourceCellIndices.at(i);
		double cellDX = intoDomain_->cellLengths[0].GetAt(cix);
		double cellDR = intoDomain_->cellLengths[1].GetAt(cix);
		double cellR = intoDomain_->localCellCenterPositions[1].GetAt(cix);
		double rInner = cellR - 0.5*cellDR;
		double rOuter = cellR + 0.5*cellDR;
		volumes.at(i) = cellDX * M_PI * (pow(rOuter, 2) - pow(rInner, 2));
	}

	double totalVolume = 0;
	for (const double& v : volumes)
		totalVolume += v;

	double averageUOutside = GetAverageFieldQuantityAroundValve(outOfDomain_->u, CURRENT_TIME_STEP, false);
	double averageVOutside = GetAverageFieldQuantityAroundValve(outOfDomain_->v, CURRENT_TIME_STEP, false);
	

	// not very pretty, but do it again. Might want to reformat this to not need re-calculation.
	for (size_t i = 0; i < sourceTermBuffer_.size(); i++)
	{
		const CellIndex& cix = sourceCellIndices.at(i);
		double cellDX = intoDomain_->cellLengths[0].GetAt(cix);
		double cellDR = intoDomain_->cellLengths[1].GetAt(cix);
		double cellR = intoDomain_->localCellCenterPositions[1].GetAt(cix);
		double rInner = cellR - 0.5*cellDR;
		double rOuter = cellR + 0.5*cellDR;

		double cellVolume = cellDX * M_PI * (pow(rOuter, 2) - pow(rInner, 2));
		double mfr_eq = totalMassFlowRate * cellVolume / totalVolume; // normalised mass flow rate

		double outerSurfaceArea = (2.0*M_PI*cellR*cellDX);
		double volumeFlowRate = mfr_eq/averageDensityOutOfDomain/outerSurfaceArea;
		sourceTermBuffer_.at(i).density = mfr_eq/cellVolume;
		sourceTermBuffer_.at(i).u = (mfr_eq*averageUOutside)/cellVolume;
		sourceTermBuffer_.at(i).v = (mfr_eq*averageVOutside)/cellVolume;
		sourceTermBuffer_.at(i).e = (gamma/(gamma-1.0)*averagePressureOutOfDomain/averageDensityOutOfDomain + 0.5*volumeFlowRate*sqrt(pow(averageUOutside, 2) + pow(averageVOutside,2)))*mfr_eq/cellVolume;
	}

	// todo: set drain terms for the other domain

	
	
}

void ReedValve::SetInitialConditions()
{
	#ifdef _DEBUG
	if (fem_.dt == 0)
		throw std::logic_error("FEM module is not initialised.");
	#endif

	// TODO: Right now creates & destroys them every time. Possible optimisation would be to cache them?
	std::vector<double> forcesOnNodes;
	std::vector<double> newDeflection; // Name based on florian's code, still need to actually figure out what it means...
	CalculateForceOnNodesFromPressure(forcesOnNodes, CURRENT_TIME_STEP);
	fem_.SolveCholeskySystem(newDeflection, forcesOnNodes);
	fem_.UpdatePositions(newDeflection);

	
}

double ReedValve::GetTipDeflection() const
{
	return fem_.nodePositionsRelativeToRoot.back().y;
}

double ReedValve::GetAverageFieldQuantityAroundValve(const FieldQuantity& fieldQuantity, const EFieldQuantityBuffer bufferName, const bool bInwards) const
{
#ifdef _DEBUG
	// First check if the given field quantity is actually on the domain that the valve is at. better safe than sorry!
	if (fieldQuantity.domain != intoDomain_ || fieldQuantity.domain != outOfDomain_)
		throw std::logic_error("Trying to get the average value of a field quantity from a valve on a domain that the valve is not connected to!");
#endif
	
	// Define a bounding box, where within all cells will be sampled.
	CellIndex boundingBoxCells[2];
	if (bInwards)
	{
		// The location in the intoDomain is saved, so that can be used directly to extrude. Extend up to the depth of the tip deflection.
		//TODO: check if normal is not flipped.
		const auto extrudedPositions = ExtrudeAlongNormal(hingePositionInDomain, holeEndPositionInDomain,GetTipDeflection());
		boundingBoxCells[0] = intoDomain_->InvertPositionToIndex(extrudedPositions.first);
		boundingBoxCells[1] = intoDomain_->InvertPositionToIndex(extrudedPositions.second);
	}
	else
	{
		// the location is in the outOfDomain. The location of the valve is not defined in this coordinate frame yet, so first transform it.
		// Alternative way to do it would be to first compute to global coordinate frame, but since it is already known the start positions are on the edge, it can be done simpler this way.
		const auto start = intoDomain_->GetLocationAlongBoundaryInAdjacentDomain(boundary_, positionAlongBoundary_);
		const Position startPos = outOfDomain_->PositionAlongBoundaryToCoordinate(start.first, start.second, 0);
		boundingBoxCells[0] = outOfDomain_->InvertPositionToIndex(startPos);
		
		const auto end = intoDomain_->GetLocationAlongBoundaryInAdjacentDomain(boundary_, holeEndPositionAlongBoundary);
		const Position endPos = outOfDomain_->PositionAlongBoundaryToCoordinate(end.first, end.second, SAMPLING_DEPTH_FOR_OUT_OF_DOMAIN);
		boundingBoxCells[1] = outOfDomain_->InvertPositionToIndex(endPos);
	}

	// Determine the total sum of all those cells
	double totalSum = 0;
	const TwoDimensionalArray& buffer = fieldQuantity.bufferMap.at(bufferName);
	for (int xIdx = boundingBoxCells[0].x; xIdx < boundingBoxCells[1].x; xIdx++)
	{
		for (int yIdx = boundingBoxCells[0].y; yIdx < boundingBoxCells[1].y; yIdx++)
		{
			totalSum+= buffer.GetAt(xIdx, yIdx);
		}
	}

	// Divide by total amount of cells to get average.
	const double sizeX = abs(boundingBoxCells[0].x - boundingBoxCells[1].x);
	const double sizeY = abs(boundingBoxCells[0].y - boundingBoxCells[1].y);
	return totalSum / (sizeX * sizeY);	
}

void ReedValve::CalculateAerodynamicDamping(std::vector<double> &forceVectorOut) //const
{
	#ifdef _DEBUG
	assert(forceVectorOut.size() == static_cast<size_t>(fem_.amountOfNodes * N_DOF_PER_NODE));
	#endif
	
	for (BeamSection& section : fem_.beamSections)
	{
		double mass = section.density * (section.b[0] + section.b[1]) * (section.h[0] + section.h[1]) * (section.length) * 0.125;

		Position posNow = fem_.GetPositionOfBeamSection(section);
		Position posPreviously = (fem_.positionsInPreviousTimeStep.at(section.leftNodeIndex) + fem_.positionsInPreviousTimeStep.at(section.rightNodeIndex)) * 0.5; // Little ugly, but not necessary to make a separate function for this.
		double leftDy = fem_.nodePositionsRelativeToRoot.at(section.leftNodeIndex).y - fem_.positionsInPreviousTimeStep.at(section.leftNodeIndex).y;
		double rightDy = fem_.nodePositionsRelativeToRoot.at(section.rightNodeIndex).y - fem_.positionsInPreviousTimeStep.at(section.rightNodeIndex).y;
		double dyDt = (posNow.y - posPreviously.y) / intoDomain_->simCase->dt;
		

		double dampingFactor; // Called epsilon in Florian (2017)
		if (dyDt >= 0)
			dampingFactor = DAMPING_C1 + DAMPING_C2 * dyDt;
		else // dy < 0
			dampingFactor = DAMPING_C3 * posNow.y;

		double dampingForce = -2 * naturalFrequency / mass * dyDt * dampingFactor;

		forceVectorOut[section.leftNodeIndex] += 0.5*dampingForce;
		forceVectorOut[section.rightNodeIndex] += 0.5*dampingForce;
	}
}

double ReedValve::FukanariReferenceArea() const
{
	if (beamProfile_ != STRAIGHT_DOUBLE_TAPERED)
	{
		throw std::logic_error("Fukanari's reference area is only confirmed for straight double tapered reed valves. Maybe the difference is small, but this is for your information. Ignore this message if you know what you're doing.");
	}
	return GetTipDeflection()*(fem_.rootWidth + sqrt(fem_.freeLength * fem_.freeLength + 0.25*pow(fem_.rootWidth - fem_.tipWidth,2)));
}

double ReedValve::DischargeCoefficient() const
{
	return fmax(0.0,fmin(1.0,0.9193*pow(1.0+0.5*GetTipDeflection()*1000.0,-0.596)));
}
