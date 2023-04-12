#include "reed_valve.h"
#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cmath>

#include "parameters.h"

#define HOLE_FACTOR 0.9

ReedValve::ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile) :
	IValve(intoDomain, boundary, positionAlongBoundary),
	FemDeformation(amountOfFreeSections, amountOfFixedNodes, beamProfile, lengthOfFreeSection, lengthOfFixedSections, intoDomain->simCase->dt, boundary),
	positionInDomain(intoDomain->PositionAlongBoundaryToCoordinate(boundary, positionAlongBoundary))
{
	holeStartPos = intoDomain->InvertPositionToIndex(intoDomain->PositionAlongBoundaryToCoordinate(boundary, positionAlongBoundary));
	holeEndPos = intoDomain->InvertPositionToIndex(intoDomain->PositionAlongBoundaryToCoordinate(boundary, positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection));
	SetSourceCellIndices(sourceCellIndices, boundary, positionAlongBoundary, lengthOfFreeSection, lengthOfFixedSections);
	//SetPressureReadingCellIndices(boundary, 2);

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

		auto cellThisSectionIsIn = intoDomain->InvertPositionToIndex(centerPos);
		double pressureGradientNormalToBeamSection = intoDomain->p.GetGradientInDirectionAndPosition(cellThisSectionIsIn,  angle);

	}
}
void ReedValve::CalculateForceOnNodes(std::vector<double>& forcesOut, const bool bAddZerosForAlignedElements) const
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
	if (bAddZerosForAlignedElements)
		forcesOut.resize(amountOfNodes * N_DOF_PER_NODE);
	else
		forcesOut.resize(amountOfNodes);
	
	
	// This iterates over the beam section elements, but we need the indices to determine the positions. Hence, up to amountOfNodes-1
	for (int nodeIdx = fixedNodes; nodeIdx < amountOfNodes - 1; nodeIdx++)
	{
		// Sample the pressures where the element is in the physical domain. To get the position where the beam section is in the total domain, the positions of the nodes that it spans between are averaged, it is converted to the reference frame of the domain (instead of the reed valve), and then added to the actual position of the valve in the domain.
		const Position beamSectionCenterPositionInDomain = positionInDomain.PlusPositionInOtherCoordinateFrame((nodePositionsRelativeToRoot[nodeIdx] + nodePositionsRelativeToRoot[nodeIdx + 1]) * 0.5); 
		const double pressureAtBeamCenter = intoDomain->p.GetInterpolatedValueAtPosition(beamSectionCenterPositionInDomain);
		const double deltaPressureWithAmbient = pressureAtBeamCenter - outOfDomain->p.At(sinkIndex);

		// We're only interested in the (locally) vertical component. hence, determine theta, the angle it makes relative to the (local) horizontal axis.
		const Position deltaPosition = nodePositionsRelativeToRoot[nodeIdx + 1] - nodePositionsRelativeToRoot[nodeIdx];
		double cosTheta = deltaPosition.x / deltaPosition.Distance(); // Just really simple pythagoras.
		
		double forceOnElement = deltaPressureWithAmbient * beamSections.at(nodeIdx).topOrBottomSurfaceArea;

		// The forces is assumed to be equally distributed over the two different nodes.
		if (bAddZerosForAlignedElements)
		{
			forcesOut[nodeIdx * N_DOF_PER_NODE] += 0.5* forceOnElement;
			forcesOut[(nodeIdx + 1) * N_DOF_PER_NODE] += 0.5* forceOnElement;
		}
		else
		{
			forcesOut[nodeIdx] += 0.5* forceOnElement;
			forcesOut[nodeIdx + 1] += 0.5* forceOnElement;
		}
		
	}
}

void ReedValve::OnRegister()
{

}

void ReedValve::SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections) const
{
	// calculate the 'starting position' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection * (1 - HOLE_FACTOR);
	auto posStart = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd);

	CellIndex sourceStartIndexOnBoundary = intoDomain->InvertPositionToIndex(posStart);		// The index (location) on the boundary where the valve starts creating a source term.
	CellIndex sourceEndIndexOnBoundary = intoDomain->InvertPositionToIndex(posEnd);		// The index (location) on the boundary where the valve stops creating a source term.

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
}
void ReedValve::GetAveragePressure() const
{
	GetAverageFieldQuantityInternal(intoDomain->p);
}
void ReedValve::Update()
{
	// TODO: Right now creates & destroys them every time. Possible optimisation would be to cache them?
	std::vector<double> forcesOnNodes;
	std::vector<double> u2Deflection; // Name based on florian's code, still need to actually figure out what it means...
	const std::vector<double> u1Deflection = {};
	CalculateForceOnNodes(forcesOnNodes, true);
	SolveCholeskySystem(u2Deflection, forcesOnNodes);
	UpdatePositions(u1Deflection, u2Deflection);
	
}
double ReedValve::GetAverageFieldQuantityInternal(const FieldQuantity &fieldQuantity) const
{
	// First check if the given field quantity is actually on the domain that the valve is at. better safe than sorry!
	if (fieldQuantity.domain != intoDomain)
		throw std::logic_error("Trying to get the average value of a field quantity from a valve on a domain that the valve is not connected to!");

	double averageValue = 0;

	// Use a bounding box, and average all the values in that.
	// todo; move bounding box depth to a higher scope
	const auto boundingBox = GetBoundingBox(4);
	const double sizeX = abs(boundingBox.first.x - boundingBox.second.x);
	const double sizeY = abs(boundingBox.first.y - boundingBox.second.y);
	for (int xIdx = boundingBox.first.x; xIdx < boundingBox.second.x; xIdx++)
	{
		for (int yIdx = boundingBox.first.y; yIdx < boundingBox.second.y; yIdx++)
		{
			averageValue+= fieldQuantity.At(xIdx, yIdx) / (sizeX * sizeY);
		}
	}

	return averageValue;	
}
std::pair<CellIndex, CellIndex> ReedValve::GetBoundingBox(const int depth) const
{
	int xOffset = 0;
	int yOffset = 0;
	
	switch (boundary)
	{
	case EBoundaryLocation::TOP:
		yOffset = -depth;
		break;
	case EBoundaryLocation::BOTTOM:
		yOffset = depth;
		break;
	case EBoundaryLocation::LEFT:
		xOffset = depth;
		break;
	case EBoundaryLocation::RIGHT:
		xOffset = -depth;
		break;
	default:
		throw std::logic_error("Invalid boundary location.");
	}
	return { holeStartPos, holeEndPos + CellIndex(xOffset, yOffset) };
}