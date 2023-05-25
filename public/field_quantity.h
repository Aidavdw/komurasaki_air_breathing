#pragma once
#include <map>
#include "2dArray.h"
#include "index2d.h"
#include "muscl.h"

// Forward Declarations
struct Position;
class Domain;

enum class EFieldQuantityBuffer
{
	CURRENT_TIME_STEP,			// The field that physically represents the value
	RUNGE_KUTTA,				// to be used for runge kutta iteration. A temporary buffer.
	NEXT_TIME_STEP,				// The value as it will be used for the next time step
	FLUX,
};

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d array.
class FieldQuantity
{
public:
	FieldQuantity() :
		domain(nullptr),
		nGhostCells(0)
	{}

	FieldQuantity(Domain* domain, const int sizeX, const int sizeY, const double initialValue = 0., const int nGhostCells=2);
	
	Domain* domain;

	int nGhostCells;

	TwoDimensionalArray currentTimeStep;			// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;			// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray nextTimeStepBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray flux;						// Gets added to nextTimeBuffer, but separately writeable for async possibility.

	// MUSCL buffers
	TwoDimensionalArray MUSCLBuffer[4]; // index with EBoundaryLocation.

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	TwoDimensionalArray& Buffer(const EFieldQuantityBuffer bufferName); // Gets read/write from buffer
	const TwoDimensionalArray& GetAtBufferConst(const EFieldQuantityBuffer bufferName) const; // Gets read/write from buffer
	
	

	// Instead of taking considering that the entire volume has one value (lumped parameter estimation), consider that the cell center has the value, and anything outside of that will be linearly interpolated between those cells.
	double GetInterpolatedValueAtPosition(const Position& atPosition, const EFieldQuantityBuffer bufferName) const;

	// Fills left/right/top/bottomFaceMUSCLBuffer variables from the given source buffer
	void PopulateMUSCLBuffers(const EFieldQuantityBuffer sourceBufferName, const double MUSCLBias, const EFluxLimiterType fluxLimiterType);
	

	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const int xIdx, const int yIdx, const EFieldQuantityBuffer buffer);

	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const CellIndex& cellIndex, const EFieldQuantityBuffer buffer)
	{
		return operator()(cellIndex.x, cellIndex.y, buffer);
	}

	inline int GetFlattenedIndex(const int xIdx, const int yIdx) const; // Helper function for getting the flattened index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].
	inline int GetFlattenedIndex(const CellIndex& cellIndex) const;
	
	double GetGradientInDirectionAndPosition(const CellIndex posIdx, const double directionAngle) const;

	double GetAverageValue(const bool bExpectUniformField) const;

	bool IsValidIndex(const CellIndex& cellIndex, const bool bAllowGhostCells) const;

private:
	int nX_ = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY_ = 0; // Amount of fields in the y-direction, not counting ghost cells.
};