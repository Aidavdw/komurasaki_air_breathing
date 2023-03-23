#include "reed_valve.h"
#include "domain.h"
#include <stdexcept>

#define HOLE_FACTOR 0.9

ReedValve::ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections) :
	FemDeformation(),
	Valve(intoDomain, boundary, positionAlongBoundary, size)
{

	SetSourceCellIndices();

	/*
	int pressureInputStartIndexOnBoundary;	// The index (location) on the boundary where measuring the average pressure to determine the FEM forces starts.
	int pressureInputEndIndexOnBoundary;	// The index (location) on the boundary where measuring the average pressure to determine the FEM forces ends.
	*/

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

void ReedValve::SetPressureReadingCellIndices(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections)
{

}
