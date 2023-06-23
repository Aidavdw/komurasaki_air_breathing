#pragma once
#include <vector>
#include "fem_deformation.h"
#include "field_quantity.h"
#include "index2d.h"
#include "IValve.h"
#include "reed_valve_geometry.h"

// forward declarations
class FieldQuantity;

// A reed valve that can bend under the loads, letting in more or less air.
class ReedValve : public IValve
{
public:
	ReedValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, const double positionAlongBoundary, const ReedValveGeometry& reedValveGeometry , const ReedValveEmpiricalParameters& reedValveEmpiricalParameters, const MaterialProperties materialProperties, const bool bMirrored, const double lengthOfFixedSections, const int amountOfFreeSections, const int amountOfFixedNodes);

	bool bMirrored;			// The reed valve always generates its geometry in the direction that goes along with the edge, preserving  aright handed coordinate system. This flag sets the geometry to be generated exactly the other way around. Note that the FEMDeformation itself is unaware of the mirroring, as it just assumes everything is in local space.
	
	double holeEndPositionAlongBoundary;	// The hole start pos and boundary are inherited from IValve, now just the end pos.
	Position hingePositionInDomain;			// The exact position in the intoDomain where the hinge is, aka where the reed valve starts.
	Position holeEndPositionInDomain;		// The exact position in the intoDomain where the hole of this reed valve ends.

	int amountOfFixedNodes;
	int amountOfFreeNodes;

	double lengthOfFreeSection;
	double lengthOfFixedSections;

protected:
	CellIndex hingePositionIndex_;		// A cached inverted position-to-cellindex of where the hinge is, aka where the reed valve starts.
	CellIndex holeEndPositionIndex_;	// A cached inverted position-to-cellindex of where the hole of this reed valve ends.
	int positionMirrorModifier_;		// As shorthand modifier that is 1 if the valve is not mirrored, or is -1 if the valve is mirrored.
	
	FemDeformation fem_;
	
	// values only saved here to be given to fem_
	ReedValveGeometry reedValveGeometry_;
	ReedValveEmpiricalParameters reedValveEmpiricalParameters_;
	MaterialProperties materialProperties_; // Right now, constant material properties for entire valve. This needs to be changed if material is anisotropic or if the material changes as a function of moving towards the tip.

public:

	// Calculates the forces in the transverse direction on all the fem sections similar to how Florian's original code did it. The optional argument adds extra zeros to get the system of equations representations.
	void CalculateForceOnNodesFromPressure(std::vector<double>& forceVectorOut, const EFieldQuantityBuffer bufferName) const;
	
	/********* Overrides from Valve interface **********/
	// Actually sets up the valve inside of the domain.
	void OnRegister() override;
	void UpdateValveState() override;
	void FillBuffer() override;

	// Reads off the field, and solves the FEM equations a first time to determine the initial deflection for the FEM members.
	void SetInitialConditions() override;
	
	double GetTipDeflection() const;

protected:
	// Used in constructor; Sets the source cell indices based on the given 
	void FillCellIndexArrayWithLine(std::vector<CellIndex>& sourceCellIndicesOut, const EFace boundary, const double positionAlongBoundary, const  double lengthOfFreeSection, const double lengthOfFixedSections) const;

	// Gets teh average field quantity in the region incscribed between the the start- and end point of the valve on the boundary and a point projected normal form this boundary. If bInwards=true, this extends into the intoDomain, if bInwards=false, this extends into the outOfDomain. 
	double GetAverageFieldQuantityAroundValve(const FieldQuantity& fieldQuantity, const EFieldQuantityBuffer bufferName, const bool bInwards = true) const;
	
	void CalculateAerodynamicDamping(std::vector<double>& forceVectorOut);

	// Calculates the approximate area that the air can pass through given a certain opening of the valve, based on Fukanari (2015)
	double FukanariReferenceArea() const;

	// Gets the discharge coeffient based on Florian (2017).
	double DischargeCoefficient() const;
	
};