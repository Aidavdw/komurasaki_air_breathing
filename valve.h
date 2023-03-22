#include "domain_enums.h"

struct Domain;


//Abstract class that represents a valve that lets air from one domain into another.
class Valve
{
	// Called to all valves upon being registered by the SimCase
	virtual void OnRegister();

	Domain* domainToSourceInto;			// The domain that the valve creates a source into.
	EBoundaryLocation boundary;			// which side of the domain the valve is attached to
	int sourceStartIndexOnBoundary;		// The index (location) on the boundary where the valve starts creating a source term.
	int sourceEndIndexOnBoundary;		// The index (location) on the boundary where the valve stops creating a source term.

	int GetSourceTermCount() const;


};