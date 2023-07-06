#pragma once
#include <stdexcept>

// Decides how this domain will be initialised
enum class EInitialisationMethod
{
	ZERO,						// Every value will be set to 0
	AMBIENT_CONDITIONS,			// Propagate the values set in SimCase, and make this domain fit those.
	FROM_INPUT_RHO_P_U_V,				// Read from input data (file or interface)
	FROM_CHAPMAN_JOUGET_SOLUTION	// Solve the Chapman Jouget solution for this domain, considering it the engine tube.
};

// Determines how the flow will interact with the boundary of the domain.
enum class EBoundaryCondition
{
	NOT_SET,
	SLIP, // slip-wall boundary condition (impermeable but not imposing zero velocity at the wall).
	NO_SLIP, // no-slip-wall boundary condition setting 0 velocity at the wall (viscous condition).
	CONNECTED, // Denotes a connection to another domain
	SUPERSONIC_INLET, // CURRENTLY UNIMPLEMENTED. Supersonic inlet condition (not able to deal with backflow). Formerly supI
	SUPERSONIC_OUTLET // CURRENT UNIMPLEMENTED. supersonic outlet condition (not able to deal with backflow). Formerly supO
};


/*  Represents a face on a rectangle; either a position, or a direction from that face, pointing towards the centre of the rectangle.
 *	Note about the definition of coordinate systems along the boundaries!
 *	Positive axes are defined to ensure that a right-handed coordinate system is maintained and that the normals of the boundary-local coordinate systems always point inwards.
 *
*
*				    TOP
*					/\
*					|
*		 + -- > --- > --- > --- +
*	L	 |			 			|		R
*	E	/\			 			/\		I
*	F	 | ->					|  ->	G
*	T	/\          /\			/\		H
*		 |          |			|		T
*		 + -- > --- > --- > --- +
*				  BOTTOM
*/
enum EFace
{
	LEFT,
	RIGHT,
	TOP,
	BOTTOM
};

// Different from a face, this shows the direction of a vector!
enum class EAxisDirection
{
	POSITIVE,
	NEGATIVE
};

// Small auxiliary function to get the name of the enum. For printing purposes.
inline std::string FaceToString(const EFace location)
{
	switch (location)
	{
	case LEFT: return {"Left"};
	case RIGHT: return {"Right"};
	case TOP: return {"Top"};
	case BOTTOM: return {"Bottom"};
	default:
		return {"Unknown location"};
	}
}

// Returns the exact opposite of the face that is given.
inline EFace Opposite(const EFace location)
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
