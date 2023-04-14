#pragma once
#include <vector>
#include "fem_deformation.h"
#include "index2d.h"
#include "IValve.h"

struct FieldQuantity;

// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public IValve
{
public:
	ReedValve(Domain* intoDomain, Domain* outOfDomain, const EBoundaryLocation boundary, const double positionAlongBoundary, const int amountOfFreeSections, const double lengthOfFreeSection, const int amountOfFixedNodes, const double lengthOfFixedSections, const EBeamProfile beamProfile, const bool bMirrored);

	bool bMirrored;			// The reed valve always generates its geometry in the direction that goes along with the edge, preserving  aright handed coordinate system. This flag sets the geometry to be generated exactly the other way around. Note that the FEMDeformation itself is unaware of the mirroring, as it just assumes everything is in local space.

	Position hingePositionInDomain;		// The exact position in the intoDomain where the hinge is, aka where the reed valve starts.
	Position holeEndPositionInDomain;	// The exact position in the intoDomain where the hole of this reed valve ends.
	
	std::vector<CellIndex> sourceCellIndices;	// The cell indices where the mass flow will be sourced into in the intoDomain.

	int amountOfFixedNodes;
	int amountOfFreeSections;

	double lengthOfFreeSection;
	double lengthOfFixedSections;

protected:
	CellIndex hingePositionIndex_;		// A cached inverted position-to-cellindex of where the hinge is, aka where the reed valve starts.
	CellIndex holeEndPositionIndex_;	// A cached inverted position-to-cellindex of where the hole of this reed valve ends.
	int positionMirrorModifier_;		// As shorthand modifier that is 1 if the valve is not mirrored, or is -1 if the valve is mirrored.
	
	EBeamProfile beamProfile_
	FemDeformation fem_;

public:

	// Currently not implemented; alternative scheme to calculate pressure, but actually normal to the fem section.
	void CalculatePressuresOnFemSections();

	// Calculates the forces in the transverse direction on all the fem sections similar to how Florian's original code did it. The optional argument adds extra zeros to get the system of equations representations.
	void CalculateForceOnNodes(std::vector<double>& forcesOut, const bool addZerosForAlignedElements) const;
	
	/********* Overrides from Valve interface **********/
	void OnRegister() override;
	void GetAveragePressure() const override;
	void Update() override;
	void SetInitialConditions() override;

protected:
	// Used in constructor; Sets the source cell indices based on the given 
	void SetSourceCellIndices(std::vector<CellIndex>& sourceCellIndicesOut, const EBoundaryLocation boundary, const double positionAlongBoundary, const  double lengthOfFreeSection, const double lengthOfFixedSections) const;
	
	double GetAverageFieldQuantityInternal(const FieldQuantity& fieldQuantity) const;

	std::pair<CellIndex, CellIndex> GetBoundingBox(const int amountOfCellsDeep=5) const;
};