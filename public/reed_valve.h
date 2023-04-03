#pragma once
#include "valve.h"
#include "fem_deformation.h"
#include <vector>

#include "index2d.h"

struct FieldQuantity;

// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public Valve, public FemDeformation
{
public:
	ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile);

	CellIndex holeStartPos; // The position (on the boundary) where the hole starts.
	CellIndex holeEndPos;	// The position (on the boundary) where the hole ends.
	
	std::vector<CellIndex> sourceCellIndices;
	//std::vector<std::pair<int, int>> pressureReadingCellIndices;


	void CalculatePressuresOnFemSections();

	void OnRegister() override;

private:
	// Used in constructor; Sets the source cell indices based on the given 
	void SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EBoundaryLocation boundary, const double positionAlongBoundary, const  double lengthOfFreeSection, const double lengthOfFixedSections) const;

	void GetAveragePressure() const override;

	void GetAverageValueNearSourceTermsInternal(const FieldQuantity& fieldQuantity) const;

	std::pair<CellIndex, CellIndex> GetBoundingBox(const int depth=5) const;
	

	// Not in use anymore; replaced by reading off the gradients directly.
	//void SetPressureReadingCellIndices(const EBoundaryLocation boundary, const int offsetFromSourceCells);
};