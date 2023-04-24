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
	Domain(const std::string& name, SimCase* simCase, const Position position, const double size[2], const int amountOfCells[2], const MeshSpacing meshSpacing[2], const EInitialisationMethod initialisationMethod);

	std::string name;
	SimCase* simCase;
	EInitialisationMethod initialisationMethod;

	Position position = {0,0};						// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};								// The total extents of the domain
	int amountOfCells[2] = {0,0};						// total amount of cells in the axis direction. This includes the ghost cells.
	std::map<EBoundaryLocation, Boundary> boundaries;	// Left, right, bottom, and up boundaries.
	MeshSpacing meshSpacing[2];							// How the grid spacing looks; first in x-direction, then in y-direction.

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


	void SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type);

	int GetTotalAmountOfCells() const;

	// Shorthand function to get the cell sizes at a certain position.
	void GetCellSizes(const CellIndex cellPos, double xSizeOut, double ySizeOut) const;

	void CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	//todo: implement output in format of CellQuantities
	//CellValues GetValuesInCell(const int xidx, const int yidx);

	// Sets all the cells to some given ambient conditions.
	void SetToAmbientConditions(const double T, const double p, const double u, const double v, const double R_ideal, const double gamma);

private:
	// Sets the domain dimensions by reference.
	void PopulateDomainDimensions();
};



void ValidateAxisInput(const int axis);