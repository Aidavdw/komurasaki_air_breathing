#pragma once
#include <ios>

// Decides how this domain will be initialised
enum EInitialisationMethod
{
	ZERO,						// Every value will be set to 0
	AMBIENT_CONDITIONS,			// Propagate the values set in SimCase, and make this domain fit those.
	FROM_INPUT_DATA,				// Read from input data (file or interface)
	FROM_CHAPMAN_JOUGET_SOLUTION	// Solve the Chapman Jouget solution for this domain, considering it the engine tube.
};

// Determines how the flow will interact with the boundary of the domain.
enum EBoundaryCondition
{
	NOT_SET,
	SLIP, // slip-wall boundary condition (impermeable but not imposing zero velocity at the wall).
	NO_SLIP, // no-slip-wall boundary condition setting 0 velocity at the wall (viscous condition).
	CONNECTED, // Denotes a connection to another domain
	SUPERSONIC_INLET, // CURRENTLY UNIMPLEMENTED. Supersonic inlet condition (not able to deal with backflow). Formerly supI
	SUPERSONIC_OUTLET // CURRENT UNIMPLEMENTED. supersonic outlet condition (not able to deal with backflow). Formerly supO
};


/*	Note about the definition of coordinate systems along the boundaries!
 *	Positive axes are defined to ensure that a right-handed coordinate system is maintained and that the normals of the boundary-local coordinate systems always point inwards.
 *
 *				    TOP
 *		 + -- > --- > --- > --- +
 *	L	 |			|			|	R
 *	E	/\			\/			\/	I
 *	F	 | ->				 <- |	G
 *	T	/\          /\			\/	H
 *		 |          |			|	T
 *		 + -- < --- < --- < --- +
 *				  BOTTOM
 *		 
 */

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