#include "reed_valve.h"
#include "domain.h"
#include "sim_case.h"
#include <stdexcept>
#include <cmath>

#define HOLE_FACTOR 0.9

ReedValve::ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile) :
	Valve(intoDomain, boundary, positionAlongBoundary),
	FemDeformation(amountOfFreeSections, amountOfFixedNodes, beamProfile, lengthOfFreeSection, lengthOfFixedSections, intoDomain->simCase.dt)
{

	SetSourceCellIndices(boundary, positionAlongBoundary, lengthOfFreeSection, lengthOfFixedSections);
	//SetPressureReadingCellIndices(boundary, 2);

}

void ReedValve::CalculatePressuresOnFemSections()
{
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
	//todo: add return based on how the code is structured.
}

void ReedValve::OnRegister()
{

}

void ReedValve::SetSourceCellIndices(const EBoundaryLocation boundary, double positionAlongBoundary, const double lengthOfFreeSection, const double lengthOfFixedSections)
{
	// calculate the 'starting position' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection * (1 - HOLE_FACTOR);
	auto posStart = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd);

	CellIndex sourceStartIndexOnBoundary = intoDomain->InvertPositionToIndex(posStart);		// The index (location) on the boundary where the valve starts creating a source term.
	CellIndex sourceEndIndexOnBoundary = intoDomain->InvertPositionToIndex(posEnd);		// The index (location) on the boundary where the valve stops creating a source term.

	// Determine all the positions between the two, and save them. Note that it can be either horizontal or vertical, so first check that
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
			sourceCellIndices.emplace_back(i, sourceStartIndexOnBoundary.y);
		}
	}
	else if (bVerticalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.y; i < sourceEndIndexOnBoundary.y; i++)
		{
			sourceCellIndices.emplace_back(sourceStartIndexOnBoundary.x, i);
		}
	}
}

/*
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

	//pressureReadingCellIndices = sourceCellIndices;

	for (auto& pos : pressureReadingCellIndices)
	{
		pos.x += xOffset;
		pos.y += yOffset;
	}
}
*/