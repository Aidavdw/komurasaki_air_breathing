#pragma once
#include "IValve.h"
#include "fem_deformation.h"
#include <vector>

#include "index2d.h"

struct FieldQuantity;

// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public IValve, public FemDeformation
{
public:
	ReedValve(Domain* intoDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile);

	Position positionInDomain;

	CellIndex holeStartPos; // The position (on the boundary) where the hole starts.
	CellIndex holeEndPos;	// The position (on the boundary) where the hole ends.
	
	std::vector<CellIndex> sourceCellIndices;
	CellIndex sinkIndex;		// The index in the outOfDomain* where the valve drains from.
	//std::vector<std::pair<int, int>> pressureReadingCellIndices;

	// Currently not implemented; alternative scheme to calculate pressure, but actually normal to the fem section.
	void CalculatePressuresOnFemSections();

	// Calculates the forces in the transverse direction on all the fem sections similar to how Florian's original code did it. The optional argument adds extra zeros to get the system of equations representations.
	void CalculateForceOnNodes(std::vector<double>& forcesOut, const bool addZerosForAlignedElements) const;

public:
	// Overrides from Valve interface
	void OnRegister() override;
	void GetAveragePressure() const override;
	void Update() override;

private:
	// Used in constructor; Sets the source cell indices based on the given 
	void SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EBoundaryLocation boundary, const double positionAlongBoundary, const  double lengthOfFreeSection, const double lengthOfFixedSections) const;



	double GetAverageFieldQuantityInternal(const FieldQuantity& fieldQuantity) const;

	std::pair<CellIndex, CellIndex> GetBoundingBox(const int amountOfCellsDeep=5) const;
	

	// Not in use anymore; replaced by reading off the gradients directly.
	//void SetPressureReadingCellIndices(const EBoundaryLocation boundary, const int offsetFromSourceCells);
};