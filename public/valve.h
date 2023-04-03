#pragma once
#include "domain_enums.h"

struct Domain;


//Abstract class that represents a valve that lets air from one domain into another.
class Valve
{
public:
	virtual ~Valve() = default;
	Valve(Domain* intoDomain, const EBoundaryLocation boundary, double positionAlongBoundary);


	EBoundaryLocation boundary;			// which side of the domain the valve is attached to
	double positionAlongBoundary;

protected:

	Domain* intoDomain;					// The domain that the valve creates a source into.

public:
	// Called to all valves upon being registered by the SimCase
	virtual void OnRegister();
	virtual void SetSourceTerms();
	virtual void GetAveragePressure() const;

private:



};