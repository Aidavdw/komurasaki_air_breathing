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
		name(name),
		XSTART(0),
		YSTART(0),
		XLENGTH(0),
		YLENGTH(0),
		GRID_RATIO_X(0),
		GRID_RATIO_Y(0),
		X_V_START(0)
	{};
	Domain(std::string& name, double XSTART, double YSTART, double XLENGTH, double YLENGTH, double GRID_RATIO_X, double GRID_RATIO_Y, double X_V_START) :
		name(name),
		XSTART(XSTART),
		YSTART(YSTART),
		XLENGTH(XLENGTH),
		YLENGTH(YLENGTH),
		GRID_RATIO_X(GRID_RATIO_X),
		GRID_RATIO_Y(GRID_RATIO_Y),
		X_V_START(X_V_START)
	{};

	std::string name;

	double XSTART;             // the x coordinate of the most bottom left point for the domains.
	double YSTART;             // the y coordinate of the most bottom left point for the domains.
	double XLENGTH;            // The length of the domains in the x-direction
	double YLENGTH;            // The length of the domains in the y-direction
	double GRID_RATIO_X;
	double GRID_RATIO_Y;
	double X_V_START;
	Boundary boundaries[4];
	SimCase* simCase = nullptr;

	void SetBoundaryType(const EBoundaryLocation location, const EBoundaryType type)
};