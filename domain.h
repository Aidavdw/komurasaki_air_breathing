#pragma once
#include <string>


class SimCase;

// Determines how the flow will interact with the boundary of the domain.
enum EBoundaryType
{
	NOT_SET,
	SLIP, // slip-wall boundary condition (impermeable but not imposing zero velocity at the wall).
	WALL, // no-slip-wall boundary condition setting 0 velocity at the wall (viscous condition).
	CONNECTED, // Denotes a connection to another domain
	SUPERSONICINLET, // Supersonic inlet condition (not able to deal with backflow). Formerly supI
	SUPERSONICOUTLET //  supersonic outlet condition (not able to deal with backflow). Formerly supO
};

enum EBoundaryLocation
{
	LEFT,
	RIGHT,
	TOP,
	BOTTOM
};

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
	Domain* domain; // The domain that this boundary is a part of
	EBoundaryType boundaryType;

	// If it's connected to another boundary too:
	Boundary* connectedBoundary;
};



// Contains information on a specific domain.
struct Domain
{
	Domain(std::string& name) :
		name(name)
	{};	

	std::string name;

	double position[2] = {0,0};			// the coordinate of the most bottom left point for the domains.
	double size[2] = {0,0};				// The total extents of the domain
	double X_V_START = 0;
	int gridResolution[2] = {0,0};		// total amount of cells in the axis direction. This includes the ghost cells.
	Boundary boundaries[4];

	void SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type);

	// Function replacing NXtot/NYtot
	int GetAmountOfCellsInAxis(const unsigned int axis) const;

	// Gets the total amount of cells in a specific axis.
	int GetCellResolutionInAxis(const int axis) const;

	int GetTotalAmountOfCells() const;
};


void ValidateAxisInput(const int axis);