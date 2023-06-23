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
	Domain(const std::string& name, SimCase* simCase, const Position& position, const std::pair<double, double> sizeArg, const std::pair<MeshSpacing, MeshSpacing> meshSpacingArg, const EInitialisationMethod initialisationMethod, const int ghostCellDepth);

	std::string name;
	SimCase* simCase;
	EInitialisationMethod initialisationMethod;

	Position position = {0,0};						// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};								// The total extents of the domain
	int amountOfCells[2] = {0,0};						// total amount of cells in the axis direction. This includes the ghost cells.
	std::map<EFace, Boundary> boundaries;	// Left, right, bottom, and up boundaries.
	MeshSpacing meshSpacing[2];							// How the grid spacing looks; first in x-direction, then in y-direction.
	int nGhost;											// How many ghost cells are generated on each side.

	FieldQuantity rho;							// Density
	FieldQuantity u;							// velocity x-component
	FieldQuantity v;							// velocity y-component
	FieldQuantity p;							// pressure
	FieldQuantity E;							// Internal energy?
	FieldQuantity T;							// Temperature
	FieldQuantity H;							// enthalpy ?

	std::vector<double> cellLengths[2];					// The length of each cell.
	std::vector<double> localCellCenterPositions[2];		// The location of the cell relative to where the domain is anchored. To get global position, add with Domain.position.

	// returns the cell indices that this position is in.
	CellIndex InvertPositionToIndex(const Position& pos) const
	{
		Position blank;
		return InvertPositionToIndex(pos, blank);
	}

	// returns the cell indices that this position is in. Also sets how far from the centre it is by reference.
	CellIndex InvertPositionToIndex(const Position pos, Position& distanceFromCenterOut) const;
	std::pair<EFace, double> GetLocationAlongBoundaryInAdjacentDomain(const EFace boundaryInThisDomain, const double positionAlongBoundaryInThisDomain) const;

	Position PositionAlongBoundaryToCoordinate(const EFace boundary, const double positionAlongBoundary, const double depth) const;

	void SetBoundaryType(const EFace location, const EBoundaryCondition type);

	int GetTotalAmountOfCells() const;
	CellIndex GetOriginIndexOfBoundary(const EFace boundary) const;
	// Gets the dimensions of the part of the ghost cells as described in the GhostOrigin reference frame.
	std::pair<int,int> GetGhostDimensions(EFace boundary);

	// Shorthand function to get the cell sizes at a certain position.
	std::pair<double, double> GetCellSizes(const CellIndex cellPos) const;
	double GetCellVolume(const CellIndex cix) const;
	double GetLengthOfSide(const EFace face) const; // Small auxiliary function that returns size[0] or size[1].

	void CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	// Sets all the cells to some given ambient conditions.
	void SetToAmbientConditions(const double temperatureSet, const double pSet, const double uSet, const double vSet);

	// calculates gamma, the specific heat ratio. Current implementation just returns a fixed value, but in reality it is dependent on species & temperature.
	double SpecificHeatRatio() const;
	double GasConstant() const;

	void UpdateGhostCells();

	// Actually do a time step. Solve fluxes, etc
	void PopulateFlowDeltaBuffer(const double dt);
	void EmptyFlowDeltaBuffer();

	void SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(const int currentRungeKuttaIter);

	bool ValidateCellIndex(const CellIndex cellIndex, const bool bAllowGhostCells) const;
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