#include "reed_valve.h"
#include "domain.h"
#include <stdexcept>

#define HOLE_FACTOR 0.9

ReedValve::ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections) :
	FemDeformation(),
	Valve(intoDomain, boundary, positionAlongBoundary, size)
{

	SetSourceCellIndices(intoDomain, boundary, positionAlongBoundary, lengthOfFreeSection, lengthOfFixedSections);
	SetPressureReadingCellIndices(boundary, 2);

}

void ReedValve::OnRegister()
{

}

void ReedValve::SetSourceCellIndices(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections)
{
	// calculate the 'starting position' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection * (1 - HOLE_FACTOR);
	auto posStart = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd);

	std::pair<int, int> sourceStartIndexOnBoundary = intoDomain->InvertPositionToIndex(posStart.first, posStart.second);		// The index (location) on the boundary where the valve starts creating a source term.
	std::pair<int, int> sourceEndIndexOnBoundary = intoDomain->InvertPositionToIndex(posEnd.first, posEnd.second);		// The index (location) on the boundary where the valve stops creating a source term.

	// Determine all the positions between the two, and save them. Note that it can be either horizontal or vertical, so first check that
	bool bHorizontalDifference = (sourceStartIndexOnBoundary.first != sourceEndIndexOnBoundary.first);
	bool bVerticalDifference = (sourceStartIndexOnBoundary.second != sourceEndIndexOnBoundary.second);

	if (bHorizontalDifference == bVerticalDifference)
		throw std::logic_error("A reed valve cannot have a source term in both directions!");

	if (bHorizontalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.first; i < sourceEndIndexOnBoundary.first; i++)
		{
			sourceCellIndices.emplace_back(i, sourceStartIndexOnBoundary.second);
		}
	}
	else if (bVerticalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.second; i < sourceEndIndexOnBoundary.second; i++)
		{
			sourceCellIndices.emplace_back(sourceStartIndexOnBoundary.first, i);
		}
	}
}

void ReedValve::SetPressureReadingCellIndices(const EBoundaryLocation boundary, const int offsetFromSourceCells)
{
	// Florian's original implementation was rather hacky and unphysical. So, instead we just choose the fields that are slightly deeper into the domain than the source cells.

	if (sourceCellIndices.size() == 0)
	{
		throw std::logic_error("Pressure reading cell index setting requires source cell indices to be set.");
	}

	int xOffset = 0;
	int yOffset = 0;
	switch (boundary)
	{
	case EBoundaryLocation::TOP:
		yOffset = -offsetFromSourceCells;
		break;
	case EBoundaryLocation::BOTTOM:
		yOffset = offsetFromSourceCells;
		break;
	case EBoundaryLocation::LEFT:
		xOffset = offsetFromSourceCells;
		break;
	case EBoundaryLocation::RIGHT:
		xOffset = -offsetFromSourceCells;
		break;
	default:
		throw std::logic_error("Invalid boundary location.");
	}

	pressureReadingCellIndices = sourceCellIndices;

	for (auto& pos : pressureReadingCellIndices)
	{
		pos.first += xOffset;
		pos.second += yOffset;
	}
}
