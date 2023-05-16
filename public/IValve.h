#pragma once
#include <vector>

#include "domain_enums.h"
#include "euler_container.h"
#include "field_quantity.h"
#include "index2d.h"

struct Domain;


//Abstract class that represents a valve that lets air from one domain into another.
class IValve
{
public:
	virtual ~IValve() = default;
	IValve(Domain* intoDomain, Domain* outOfDomain, const EFace boundary, double positionAlongBoundary);


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

	void CalculateFlow();					// Safely calls FillBuffer. Proxy function.
	
	virtual void OnRegister();				// Called to all valves upon being registered by the SimCase
	virtual void UpdateValveState();		// Samples the environment, calculates forces etc, and updates how open a valve is.
	virtual void SetInitialConditions();	

	virtual void AddBufferTermsToSourceCells(const EFieldQuantityBuffer bufferToWriteInto);	// Sets the source terms in sourceTermBuffer_ as calculated with CalculateFlowFromValve(), and puts them in the sourceCellsIndices in the intoDomain_. Excepts if it cannot do this.
	virtual void EmptyBuffer();				// Sets the sourceTermBuffer and sinkTermBuffer back to zeros.

protected:
	virtual void FillBuffer(); // Don't call manually, this is supposed to be called in proxy by CalculateFlow(). Calculates the flow from the valve, and fills it into the sourceTermBuffer_ variable.

};