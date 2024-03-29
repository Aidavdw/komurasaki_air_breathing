#pragma once
#include <vector>

#include "domain_enums.h"
#include "euler_container.h"
#include "field_quantity.h"
#include "index2d.h"

// forward declarations
class Domain;


//Abstract class that represents a valve that lets air from one domain into another. This is how the simulation interacts with valves, so any valve must be derived from this class.
class IValve
{
public:
	virtual ~IValve();
	IValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, double positionAlongBoundary);

	Position pos;
	std::string label;


protected:
	EFace boundary_;										// which side of the domain the valve is attached to
	double positionAlongBoundary_;
	Domain* outOfDomain_;									// the domain that the valve sinks out of
	Domain* intoDomain_;									// The domain that the valve creates a source into.
	std::vector<CellIndex> sourceCellsIndices_;				// The cells in the intoDomain that this valve will add an accumulation term to.
	std::vector<EulerContinuity> sourceTermBuffer_;			// Stores how much will be sourced into the sourceCellIndices.

	std::vector<CellIndex> sinkCellsIndices_;				// The cells in the outOfDomain that this valve will subtract the accumulation term from.
	std::vector<EulerContinuity> sinkTermBuffer_;			// Stores how much will be sourced into the sourceCellIndices.


public:

	void CalculateFlow();					// Safely calls FillBuffer(). Proxy function.
	std::string ToString() const;
	void AddCachedTermsToDomainConservationEquations() const;	// Uses the values cached into sourceTermBuffer_ as calculated with CalculateFlowFromValve(), and puts them in the sourceCellsIndices in the intoDomain_. Excepts if it cannot do this.
	void EmptyBuffer();				// Sets the sourceTermBuffer and sinkTermBuffer back to zeros.


/********** TO IMPLEMENT INTERFACE **********/
public:
	virtual void OnRegister();				// Called to all valves upon being registered by the SimCase
	virtual void UpdateValveState();		// Samples the environment, calculates forces etc, and updates how open a valve is.
	virtual void SetInitialConditions();	// Called before first 

protected:
	virtual void FillSourceTermBuffer(); // Don't call manually, this is supposed to be called in proxy by CalculateFlow(). Calculates the flow from the valve, and fills it into the sourceTermBuffer_ variable.

};