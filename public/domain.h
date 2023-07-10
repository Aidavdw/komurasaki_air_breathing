#pragma once
#include <string>
#include "field_quantity.h"
#include "mesh_spacing.h"
#include "pos2d.h"
#include "index2d.h"
#include <map>

#include "boundary.h"
#include "domain_enums.h"


class SimCase;


// Represents a volume containing fluid. The properties of the fluid at any position are described by its FieldQuantities
class Domain
{
public:
	Domain(const std::string& name, SimCase* simCase, const Position& position, const std::pair<double, double> sizeArg, const std::pair<MeshSpacingSolution, MeshSpacingSolution> meshSpacingArg, const EInitialisationMethod initialisationMethod, const int ghostCellDepth);

	std::string name;
	SimCase* simCase;
	EInitialisationMethod initialisationMethod;

	Position position = {0,0};						// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};								// The total extents of the domain
	int amountOfCells[2] = {0,0};						// total amount of cells in the axis direction. This includes the ghost cells.
	std::map<EFace, Boundary> boundaries;	// Left, right, bottom, and up boundaries.
	MeshSpacingSolution meshSpacing[2];							// How the grid spacing looks; first in x-direction, then in y-direction.
	int nGhost;											// How many ghost cells are generated on each side.

	FieldQuantity rho;							// Density
	FieldQuantity u;							// velocity x-component
	FieldQuantity v;							// velocity y-component
	FieldQuantity p;							// pressure
	FieldQuantity E;							// Internal energy?
	FieldQuantity T;							// Temperature
	FieldQuantity H;							// enthalpy ?

	TwoDimensionalArray eulerConservationTerms[4];		// To be used to calculate the values in the next time steps, contains the mass, momentum and energy fluxes for all terms. Written indirectly so that IValve can easily add here.

	std::vector<double> cellLengths[2];					// The length of each cell.
	std::vector<double> localCellCenterPositions[2];		// The location of the cell relative to where the domain is anchored. To get global position, add with Domain.position.

	CellIndex InvertPositionToIndex(const Position pos, Position& distanceFromCenterOut) const; 	// returns the cell indices that this position is in. Also sets how far from the centre it is by reference.
	CellIndex InvertPositionToIndex(const Position& pos) const // Proxy where the distance from centre out is discarded.
	{
		Position blank(0,0);
		return InvertPositionToIndex(pos, blank);
	}
	

	std::pair<EFace, double> GetLocationAlongBoundaryInAdjacentDomain(const EFace boundaryInThisDomain, const double positionAlongBoundaryInThisDomain) const; // Expresses the position in terms of its complementary boundary.
	Position PositionAlongBoundaryToCoordinate(const EFace boundary, const double positionAlongBoundary, const double depth) const; // Transforms a given distance out from the datum on any boundary to a origin-relative coordinate.

	void SetBoundaryType(const EFace location, const EBoundaryCondition type);

	int GetTotalAmountOfCells() const;
	CellIndex GetOriginIndexOfBoundary(const EFace boundary) const; // Returns the datum of where the 'zero' cell is for a certain boundary.
	std::pair<int,int> GetGhostDimensions(EFace boundary); // Gets the dimensions of the part of the ghost cells as described in the GhostOrigin reference frame.
	std::pair<double, double> GetCellSizes(const CellIndex cellPos) const; // Shorthand function to get the cell sizes at a certain position.
	double GetCellVolume(const CellIndex cix) const;
	double GetLengthOfSide(const EFace face) const; // Small auxiliary function that returns size[0] or size[1].

	void CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to); // Copies all the entries from one buffer into another. This leaves the original intact.
	void SetToAmbientConditions(const double temperatureSet, const double pSet, const double uSet, const double vSet); // Sets all the cells to some given ambient conditions.

	// calculates gamma, the specific heat ratio. Current implementation just returns a fixed value, but in reality it is dependent on species & temperature.
	double SpecificHeatRatio() const;
	double GasConstant() const;

	void UpdateGhostCells(); // Sets the values of the ghost cells based on the type of boundary condition this has.
	void CacheEulerConservationTerms(const double dt); // Calculate the fluxes between cells, and populate the euler conservation matrices.
	void EmptyFlowDeltaBuffer();
	void SetNextTimeStepValuesBasedOnCachedEulerContinuities(const int currentRungeKuttaIter);

	bool ValidateCellIndex(const CellIndex cellIndex, const bool bAllowGhostCells) const; // Checks whether this cell index is actually valid and/or inside of this domain. Mostly used for debugging.
	bool AreAllBoundariesSet() const; // Checks if all boundaries are set up
	

private:
	// Caches the domain dimensions, and saves them in cellLengths and localCellCenterPositions.
	void CacheCellSizes();

	void PopulateSlipConditionGhostCells(const EFace boundary);
	void PopulateNoSlipConditionGhostCells(const EFace boundary);
	void PopulateConnectedGhostCells(const EFace boundary);
	//void PopulateSupersonicInletGhostCells(const EFace boundary);
	//void PopulateSupersonicOutletGhostCells(const EFace boundary);
};



void ValidateAxisInput(const int axis);