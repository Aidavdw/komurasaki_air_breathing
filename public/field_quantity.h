#pragma once
#include <map>
#include "2dArray.h"
#include "index2d.h"
#include "muscl.h"

struct Position;
struct Domain;

enum EFieldQuantityBuffer
{
	CURRENT_TIME_STEP,			// The field that physically represents the value
	RUNGE_KUTTA,				// to be used for runge kutta iteration. A temporary buffer.
	NEXT_TIME_STEP,				// not yet sure what this means but oh well!
	DELTA_FLOW,
	DELTA_VALVE
};

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d array.
struct FieldQuantity
{
	FieldQuantity() :
		domain(nullptr),
		nGhostCells(0)
	{}

	FieldQuantity(Domain* domain, const int sizeX, const int sizeY, const double initialValue = 0., const int nGhostCells=2);

	// TODO: check if this actually references the desired buffers, and does not make a new copy inline!
	std::map<EFieldQuantityBuffer, TwoDimensionalArray&> bufferMap;
	Domain* domain;

	int nGhostCells;

	TwoDimensionalArray currentTimeStep;			// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;			// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray nextTimeStepBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray deltaDueToFlow;				// Gets added to nextTimeBuffer, but separately writeable for async possibility.
	TwoDimensionalArray deltaDueToValve;			// Gets added to nextTimeBuffer, but separately writeable for async possibility.

	// MUSCL buffers
	TwoDimensionalArray MUSCLBuffer[4]; // index with EBoundaryLocation.

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo);

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	// Instead of taking considering that the entire volume has one value (lumped parameter estimation), consider that the cell center has the value, and anything outside of that will be linearly interpolated between those cells.
	double GetInterpolatedValueAtPosition(const Position& atPosition, const EFieldQuantityBuffer bufferName) const;

	// Fills left/right/top/bottomFaceMUSCLBuffer variables from the given source buffer
	void PopulateMUSCLBuffers(const EFieldQuantityBuffer sourceBuffer, const double MUSCLBias, const EFluxLimiterType fluxLimiterType);
	

	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const int xIdx, const int yIdx, const EFieldQuantityBuffer buffer)
	{
		auto& buf = bufferMap.at(buffer);
		return buf[GetFlattenedIndex(xIdx, yIdx)];
	}

	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const CellIndex& cellIndex, const EFieldQuantityBuffer buffer)
	{
		auto& buf = bufferMap.at(buffer);
		return buf[GetFlattenedIndex(cellIndex.x, cellIndex.y)];
	}

	/* Should defer all this just to the buffer names.
	// Shorthand overloaded accessor for currentTimeStep buffer.
	inline double& operator () (const int xIdx, const int yIdx)
	{
		return currentTimeStep[At(xIdx, yIdx)];
	}

	// Shorthand overloaded accessor for currentTimeStep buffer.
	inline double& operator () (const CellIndex& cellIndex)
	{
		return currentTimeStep[At(cellIndex.x, cellIndex.y)];
	}
	

	inline double GetAt(const CellIndex& cellIndex) const;
	*/

	inline int GetFlattenedIndex(const int xIdx, const int yIdx) const; // Helper function for getting the flattened index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].
	inline int GetFlattenedIndex(const CellIndex& cellIndex) const;
	
	double GetGradientInDirectionAndPosition(const CellIndex posIdx, const double directionAngle) const;

	double GetAverageValue(const bool bExpectUniformField) const;

private:
	int nX_ = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY_ = 0; // Amount of fields in the y-direction, not counting ghost cells.
};