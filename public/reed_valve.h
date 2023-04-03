#pragma once
#include "valve.h"
#include "fem_deformation.h"
#include <vector>
#include "index2d.h"


// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public Valve, public FemDeformation
{
public:
	ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile);

	// TODO: Replace all the weird int[2] and double[2] with IntCoordinate and DoubleCoordinate.
	std::vector<CellIndex> sourceCellIndices;
	//std::vector<std::pair<int, int>> pressureReadingCellIndices;


	void CalculatePressuresOnFemSections();

	void OnRegister() override;

private:
	void SetSourceCellIndices(const EBoundaryLocation boundary, const double positionAlongBoundary, const  double lengthOfFreeSection, const double lengthOfFixedSections);

	// Not in use anymore; replaced by reading off the gradients directly.
	//void SetPressureReadingCellIndices(const EBoundaryLocation boundary, const int offsetFromSourceCells);
};