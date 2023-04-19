#include "reed_valve.h"

#include <cassert>

#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cmath>

#include "parameters.h"

#define HOLE_FACTOR 0.9

#define DAMPING_C1 5.0E-8 //5.0E-7  
#define DAMPING_C2 0.0E-8 //2.0E-8 // Why on earth is this 0??
#define DAMPING_C3 0.0007

ReedValve::ReedValve(Domain* intoDomain, Domain* outOfDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile, const bool bMirrored) :
	IValve(intoDomain,outOfDomain,boundary,positionAlongBoundary),
	bMirrored(bMirrored),
	amountOfFixedNodes(amountOfFixedNodes),
	amountOfFreeNodes(amountOfFreeSections),
	lengthOfFreeSection(lengthOfFreeSection),
	lengthOfFixedSections(lengthOfFixedSections),
	beamProfile_(beamProfile)
{
	positionMirrorModifier_ = bMirrored ? -1 : 1;
	hingePositionInDomain = intoDomain->PositionAlongBoundaryToCoordinate(boundary,positionAlongBoundary,0);
	holeEndPositionInDomain = intoDomain->PositionAlongBoundaryToCoordinate(boundary,positionAlongBoundary+(lengthOfFixedSections+lengthOfFreeSection)*positionMirrorModifier_,0);

	hingePositionIndex_ = intoDomain->InvertPositionToIndex(hingePositionInDomain);
	holeEndPositionIndex_ = intoDomain->InvertPositionToIndex(holeEndPositionInDomain);

	fem_ = FemDeformation(amountOfFreeSections,amountOfFixedNodes,beamProfile,lengthOfFreeSection,lengthOfFixedSections,intoDomain->simCase->dt,boundary);
}

void ReedValve::CalculatePressuresOnFemSections()
{
	throw std::logic_error("This function is not implemented yet.");
	for (int beamIdx = 0; beamIdx < beamSections.size(); beamIdx++)
	{
		// the 'left' position
		const double& pos1X = nodePositionsRelativeToRoot[beamIdx][0];
		const double& pos1Y = nodePositionsRelativeToRoot[beamIdx][1];

		// the 'right' position
		const double& pos2X = nodePositionsRelativeToRoot[beamIdx][0];
		const double& pos2Y = nodePositionsRelativeToRoot[beamIdx][1];

		// Hence, the centre, and the angle this section is at
		Position centerPos = { (pos1X + pos2X) * 0.5, (pos1Y + pos2Y) * 0.5 };

		double angle = atan2(pos2Y - pos1Y, pos2X - pos1X);

		auto cellThisSectionIsIn = intoDomain_->InvertPositionToIndex(centerPos);
		double pressureGradientNormalToBeamSection = intoDomain_->p.GetGradientInDirectionAndPosition(cellThisSectionIsIn,  angle);

	}
}
void ReedValve::CalculateForceOnNodes(std::vector<double>& forceVectorOut) const
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

	forceVectorOut.resize(fem_.amountOfNodes * N_DOF_PER_NODE);
	
	// This iterates over the beam section elements, but we need the indices to determine the positions. Hence, up to amountOfNodes-1
	for (int nodeIdx = fem_.fixedNodes; nodeIdx < fem_.amountOfNodes - 1; nodeIdx++)
	{
		// Sample the pressures where the element is in the physical domain. To get the position where the beam section is in the total domain, the positions of the nodes that it spans between are averaged, it is converted to the reference frame of the domain (instead of the reed valve), and then added to the actual position of the valve in the domain.
		const Position beamSectionCenterPositionLocal = (fem_.nodePositionsRelativeToRoot[nodeIdx] + fem_.nodePositionsRelativeToRoot[nodeIdx + 1]) * 0.5;
		const Position beamSectionCenterPositionInDomain = TransformToOtherCoordinateSystem(beamSectionCenterPositionLocal, hingePositionInDomain, {0,0});
		const double pressureAtBeamCenter = intoDomain_->p.GetInterpolatedValueAtPosition(beamSectionCenterPositionInDomain);

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

		double pressureAtSink =  outOfDomain_->p.GetInterpolatedValueAtPosition(ambientSamplePositionLocal);
		
		const double deltaPressureWithAmbient = pressureAtBeamCenter - pressureAtSink;

		// We're only interested in the (locally) vertical component. hence, determine theta, the angle it makes relative to the (local) horizontal axis.
		const Position deltaPosition = fem_.nodePositionsRelativeToRoot[nodeIdx + 1] - fem_.nodePositionsRelativeToRoot[nodeIdx];
		double cosTheta = deltaPosition.x / deltaPosition.Distance(); // Just really simple pythagoras.
		
		double forceOnElement = deltaPressureWithAmbient * fem_.beamSections.at(nodeIdx).topOrBottomSurfaceArea * cosTheta;

		// The forces is assumed to be equally distributed over the two different nodes.
		forceVectorOut[nodeIdx * N_DOF_PER_NODE] += 0.5* forceOnElement;
		forceVectorOut[(nodeIdx + 1) * N_DOF_PER_NODE] += 0.5* forceOnElement;
		
	}
}

void ReedValve::OnRegister()
{
	fem_ = FemDeformation(amountOfFreeNodes, amountOfFixedNodes, beamProfile_, lengthOfFreeSection, lengthOfFixedSections, intoDomain_->simCase->dt, boundary_);
	SetSourceCellIndices(sourceCellIndices, boundary_, positionAlongBoundary_, lengthOfFreeSection, lengthOfFixedSections);
	
}

void ReedValve::SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections) const
{
	// calculate the 'starting position' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection * (1 - HOLE_FACTOR);
	auto posStart = intoDomain_->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain_->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd);

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
void ReedValve::GetAveragePressure() const
{
	GetAverageFieldQuantityInternal(intoDomain_->p);
}

void ReedValve::Update()
{
	#ifdef _DEBUG
	if (fem_.dt == 0)
		throw std::logic_error("FEM module is not initialised.");
	#endif
	
	// TODO: Right now creates & destroys them every time. Possible optimisation would be to cache them?
	std::vector<double> forcesOnNodes;
	std::vector<double> newDeflection;
	CalculateForceOnNodes(forcesOnNodes);
	CalculateAerodynamicDamping(forcesOnNodes);
	fem_.CalculateNewDeflections(newDeflection, forcesOnNodes);
	fem_.UpdatePositions(newDeflection);
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
	CalculateForceOnNodes(forcesOnNodes);
	fem_.SolveCholeskySystem(newDeflection, forcesOnNodes);
	fem_.UpdatePositions(newDeflection);

	
}

double ReedValve::GetAverageFieldQuantityInternal(const FieldQuantity &fieldQuantity) const
{
	#ifdef _DEBUG
	// First check if the given field quantity is actually on the domain that the valve is at. better safe than sorry!
	if (fieldQuantity.domain != intoDomain_)
		throw std::logic_error("Trying to get the average value of a field quantity from a valve on a domain that the valve is not connected to!");
	#endif

	double totalSum = 0;

	// Use a bounding box, and average all the values in that.
	// todo; move bounding box depth to a higher scope
	const auto boundingBox = GetBoundingBox(4);
	for (int xIdx = boundingBox.first.x; xIdx < boundingBox.second.x; xIdx++)
	{
		for (int yIdx = boundingBox.first.y; yIdx < boundingBox.second.y; yIdx++)
		{
			totalSum+= fieldQuantity.At(xIdx, yIdx);
		}
	}

	const double sizeX = abs(boundingBox.first.x - boundingBox.second.x);
	const double sizeY = abs(boundingBox.first.y - boundingBox.second.y);
	return totalSum / (sizeX * sizeY);	
}
std::pair<CellIndex, CellIndex> ReedValve::GetBoundingBox(const int amountOfCellsDeep) const
{
	int xOffset = 0;
	int yOffset = 0;
	
	switch (boundary_)
	{
	case EBoundaryLocation::TOP:
		yOffset = -amountOfCellsDeep;
		break;
	case EBoundaryLocation::BOTTOM:
		yOffset = amountOfCellsDeep;
		break;
	case EBoundaryLocation::LEFT:
		xOffset = amountOfCellsDeep;
		break;
	case EBoundaryLocation::RIGHT:
		xOffset = -amountOfCellsDeep;
		break;
	default:
		throw std::logic_error("Invalid boundary location.");
	}
	return { hingePositionIndex_, holeEndPositionIndex_ + CellIndex(xOffset, yOffset) };
}
void ReedValve::CalculateAerodynamicDamping(std::vector<double> &forceVectorOut) //const
{
	#ifdef _DEBUG
	assert(forceVectorOut.size() == fem_.amountOfNodes * N_DOF_PER_NODE);
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