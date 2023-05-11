#pragma once
#include <string>
#include "field_quantity.h"
#include "mesh_spacing.h"
#include "pos2d.h"
#include "index2d.h"
#include <map>

#include "boundary.h"
#include "domain_enums.h"


struct SimCase;


// Contains information on a specific domain.
struct Domain
{
	Domain(const std::string& name, SimCase* simCase, const Position& position, const double sizeArg[2], const int amountOfCellsArg[2], const MeshSpacing meshSpacingArg[2], const EInitialisationMethod initialisationMethod, const int ghostCellDepth);

	std::string name;
	SimCase* simCase;
	EInitialisationMethod initialisationMethod;

	Position position = {0,0};						// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};								// The total extents of the domain
	int amountOfCells[2] = {0,0};						// total amount of cells in the axis direction. This includes the ghost cells.
	std::map<EBoundaryLocation, Boundary> boundaries;	// Left, right, bottom, and up boundaries.
	MeshSpacing meshSpacing[2];							// How the grid spacing looks; first in x-direction, then in y-direction.
	int nGhost;											// How many ghost cells are generated on each side.

	FieldQuantity rho;							// Density
	FieldQuantity u;							// velocity x-component
	FieldQuantity v;							// velocity y-component
	FieldQuantity p;							// pressure
	FieldQuantity E;							// Internal energy?
	FieldQuantity T;							// Temperature
	FieldQuantity H;							// enthalpy ?

	FieldQuantity cellLengths[2];					// The length of each cell.
	FieldQuantity localCellCenterPositions[2];		// The location of the cell relative to where the domain is anchored. To get global position, add with Domain.position.

	// returns the cell indices that this position is in.
	CellIndex InvertPositionToIndex(const Position pos) const;
	CellIndex InvertPositionToIndex(const Position pos, Position& distanceFromCenter) const;
	std::pair<EBoundaryLocation, double> GetLocationAlongBoundaryInAdjacentDomain(const EBoundaryLocation boundaryInThisDomain, const double positionAlongBoundaryInThisDomain) const;

	Position PositionAlongBoundaryToCoordinate(const EBoundaryLocation boundary, const double positionAlongBoundary, const double depth) const;


	void SetBoundaryType(const EBoundaryLocation location, const EBoundaryCondition type);
	void ConnectBoundary(const EBoundaryLocation location, Domain* otherDomain);

	int GetTotalAmountOfCells() const;
	CellIndex GetOriginIndexOfBoundary(const EBoundaryLocation boundary) const;
	// Gets the dimensions of the part of the ghost cells as described in the GhostOrigin reference frame.
	std::pair<int,int> GetGhostDimensions(EBoundaryLocation boundary);

	// Shorthand function to get the cell sizes at a certain position.
	std::pair<double, double> GetCellSizes(const CellIndex cellPos) const;

	void CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	//todo: implement output in format of CellQuantities
	//CellValues GetValuesInCell(const int xIdx, const int yIdx);

	// Sets all the cells to some given ambient conditions.
	void SetToAmbientConditions(const double temperatureSet, const double pSet, const double uSet, const double vSet, const double R_ideal, const double gamma);

	// calculates gamma, the specific heat ratio. Current implementation just returns a fixed value, but in reality it is dependent on species & temperature.
	double SpecificHeatRatio() const;
	double GasConstant() const;

	void UpdateGhostCells();

	// Actually do a time step. Solve fluxes, etc
	void PopulateFlowDeltaBuffer(const double dt, const int currentRungeKuttaIter);

	void SetNextTimeStepValuesBasedOnRungeKuttaAndDeltaBuffers(const int currentRungeKuttaIter);
	

private:
	// Caches the domain dimensions, and asves them in cellLengths and localCellCenterPositions.
	void PopulateDomainDimensions();

	void PopulateSlipConditionGhostCells(const EBoundaryLocation boundary);
	void PopulateNoSlipConditionGhostCells(const EBoundaryLocation boundary);
	void PopulateConnectedGhostCells(const EBoundaryLocation boundary);
	void PopulateSupersonicInletGhostCells(const EBoundaryLocation boundary);
	void PopulateSupersonicOutletGhostCells(const EBoundaryLocation boundary);

	bool ValidateCellIndex(const CellIndex cellIndex, const bool bAllowGhostCells) const;
};



void ValidateAxisInput(const int axis);