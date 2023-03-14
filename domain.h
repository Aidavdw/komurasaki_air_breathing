#pragma once
#include <string>
#include "field_quantity.h"
#include "mesh_spacing.h"

#include "domain_enums.h"


class SimCase;

// Describes a boundary that a domain has four of (considering it is 2d, and all domains are rectangular
struct Boundary
{
	Boundary() :
		boundaryType(EBoundaryType::NOT_SET),
		connectedBoundary(nullptr)
	{};

	Boundary(EBoundaryType boundaryType) :
		boundaryType(boundaryType),
		connectedBoundary(nullptr)
	{};

	// This boundary itself:
	//Domain* domain; // The domain that this boundary is a part of
	EBoundaryType boundaryType;

	// If it's connected to another boundary too:
	Boundary* connectedBoundary;
};



// Contains information on a specific domain.
struct Domain
{
	Domain(std::string& name, const double position[2], const double size[2], const int amountOfCells[2], const MeshSpacing meshSpacing[2], const EInitialisationMethod initialisationMethod);

	std::string name;
	EInitialisationMethod initialisationMethod;

	double position[2] = {0,0};					// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};						// The total extents of the domain
	double X_V_START = 0;
	int amountOfCells[2] = {0,0};				// total amount of cells in the axis direction. This includes the ghost cells.
	Boundary boundaries[4];						// Left, right, bottom, and up boundaries. Access using EBoundaryLocation struct.
	MeshSpacing meshSpacing[2];								

	FieldQuantity rho; // Density
	FieldQuantity u; // velocity x-component
	FieldQuantity v; // velocity y-component
	FieldQuantity p; // pressure
	FieldQuantity E; // Internal energy?
	FieldQuantity T; // Temperature
	FieldQuantity H; // enthalpy ?

	void SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type);

	// Gets the total amount of cells in a specific axis.
	int GetCellResolutionInAxis(const int axis) const;

	int GetTotalAmountOfCells() const;

	void CopyFieldQuantitiesToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	//todo: implement output in format of CellQuantities
	//CellValues GetValuesInCell(const int xidx, const int yidx);

	// Sets all the cells to some given ambient conditions.
	void SetToAmbientConditions(const double T, const double p, const double u, const double v, const double R_ideal, const double gamma);
};



void ValidateAxisInput(const int axis);