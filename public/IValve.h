#pragma once
#include "domain_enums.h"

struct Domain;


//Abstract class that represents a valve that lets air from one domain into another.
class IValve
{
public:
	virtual ~IValve() = default;
	IValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, double positionAlongBoundary);




protected:
	EFace boundary_;			// which side of the domain the valve is attached to
	double positionAlongBoundary_;
	Domain* outOfDomain_;				// the domain that the valve sinks out of
	Domain* intoDomain_;					// The domain that the valve creates a source into.

public:
	// Called to all valves upon being registered by the SimCase
	virtual void OnRegister();
	virtual void Update();

	// Gets the mass flow rate going through the valve in its current condition. Positive values ar from the outOfDomain_ into the intoDomain_
	virtual double GetMassFlowRate() const;
	virtual void SetInitialConditions();
	virtual void PopulateValveDeltaBuffer();

private:



};