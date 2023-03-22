#include "reed_valve.h"
#include "domain.h"
#include <stdexcept>

#define HOLE_FACTOR 0.9

ReedValve::ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections) :
	FemDeformation(),
	Valve(intoDomain, boundary, positionAlongBoundary, size)
{
	// calculate the 'starting positiong' based on the position along the boundary as provided, offsetting with the hole size etc.
	double posAlongBoundaryStart = positionAlongBoundary + lengthOfFixedSections +  lengthOfFreeSection*(1-HOLE_FACTOR);
	auto posStart = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryStart);
	double posAlongBoundaryEnd = positionAlongBoundary + lengthOfFixedSections + lengthOfFreeSection;
	auto posEnd = intoDomain->PositionAlongBoundaryToCoordinate(boundary, posAlongBoundaryEnd);

	std::pair<int,int> sourceStartIndexOnBoundary = intoDomain->InvertPositionToIndex(posStart.first, posStart.second);		// The index (location) on the boundary where the valve starts creating a source term.
	std::pair<int, int> sourceEndIndexOnBoundary = intoDomain->InvertPositionToIndex(posEnd.first, posEnd.second);		// The index (location) on the boundary where the valve stops creating a source term.

	// Determine all the positions between the two, and save them. Note that it can be either horizontal or vertical, so first check that
	bool bHorizontalDifference = (sourceStartIndexOnBoundary.first != sourceEndIndexOnBoundary.first);
	bool bVerticalDifference = (sourceStartIndexOnBoundary.second != sourceEndIndexOnBoundary.second);

	if (bHorizontalDifference == bVerticalDifference)
		throw std::logic_error("A reed valve cannot have a source term in both directions!");

	std::vector<double[2]> sourceCells;

	if (bHorizontalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.first; i < sourceEndIndexOnBoundary.first; i++)
		{
			sourceCells.emplace_back( i, sourceStartIndexOnBoundary.second );
		}
	}
	else if (bVerticalDifference)
	{
		for (int i = sourceStartIndexOnBoundary.second; i < sourceEndIndexOnBoundary.second; i++)
		{
			sourceCells.emplace_back(sourceStartIndexOnBoundary.first, i);
		}
	}


	/*
	int pressureInputStartIndexOnBoundary;	// The index (location) on the boundary where measuring the average pressure to determine the FEM forces starts.
	int pressureInputEndIndexOnBoundary;	// The index (location) on the boundary where measuring the average pressure to determine the FEM forces ends.

	// Sets the x-index where the source term cells due to the presence of a valve start. 
	while ((mfr_index_inf[k] < NXtot[dom_low] - NGHOST) && (x[dom_low][mfr_index_inf[k]][NGHOST] < X_V_START[k] + L_FIX + L_HOLE * (1.0 - HOLE_FACTOR)))
	{
		mfr_index_inf[k]++;
	}
	// Sets the x-index where the source term cells due to the presence of a valve end.
	while (mfr_index_sup[k] < NXtot[dom_low] - NGHOST && x[dom_low][mfr_index_sup[k] + 1][NGHOST] < X_V_START[k] + L_FIX + L_HOLE)
	{
		mfr_index_sup[k]++;
	}
	// the amount of source cells that this valve has, between the start- and end.
	mfr_n[k] = 1 + mfr_index_sup[k] - mfr_index_inf[k];
	*/

}

void ReedValve::OnRegister()
{

}
