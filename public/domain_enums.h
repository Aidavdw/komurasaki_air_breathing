#pragma once
#include <ios>

// Decides how this domain will be initialised
enum EInitialisationMethod
{
	ZERO,						// Every value will be set to 0
	AMBIENTCONDITIONS,			// Propogate the values set in SimCase, and make this domain fit those.
	FROMINPUTDATA,				// Read from input data (file or interface)
	FROMCHAPMANJOUGETSOLUTION	// Solve the Chapman Jouget solution for this domain, considering it the engine tube.
};

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

inline EBoundaryLocation Opposite(const EBoundaryLocation location)
{
	switch (location)
	{
	case TOP:
		return BOTTOM;
	case BOTTOM:
		return TOP;
	case LEFT:
		return RIGHT;
	case RIGHT:
		return LEFT;
	default:
			throw std::logic_error("Handling the Opposite() for this side is not implemented.");
	}
}