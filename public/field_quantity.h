#pragma once
#include "2dArray.h"
#include "index2d.h"

// Forward Declarations
struct Position;
class Domain;

enum class EFieldQuantityBuffer
{
	CURRENT_TIME_STEP,			// The field that physically represents the value
	RUNGE_KUTTA,				// to be used for runge kutta iteration. A temporary buffer.
	NEXT_TIME_STEP,				// The value as it will be used for the next time step
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
	
	Domain* domain;	// Domain that this field quantity is a part of.
	int nGhostCells;	// The amount of ghost cells this field quantity has.

	TwoDimensionalArray currentTimeStep;			// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;			// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray nextTimeStepBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!
	
	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to); // Copies all the values of one buffer to the other. Expensive operation!
	TwoDimensionalArray& Buffer(const EFieldQuantityBuffer bufferName); // Gets read-only from buffer
	const TwoDimensionalArray& GetAtBufferConst(const EFieldQuantityBuffer bufferName) const; // Gets read/write from buffer
	
	double GetInterpolatedValueAtPosition(const Position& atPosition, const EFieldQuantityBuffer bufferName) const; // Instead of taking considering that the entire volume has one value (lumped parameter estimation), consider that the cell center has the value, and anything outside of that will be linearly interpolated between those cells.
	double GetGradientInDirectionAndPosition(const CellIndex posIdx, const double directionAngle) const; // returns a directional derivative. //todo: note that this is not weighted for cell size. Implement in Domain layer instead.
	double GetAverageValue(const bool bExpectUniformField) const;	// Returns the average value of the entire domain. if bExpectUniformField is true, then it just samples a few points to save computation time.
	
	bool IsValidIndex(const CellIndex& cellIndex, const bool bAllowGhostCells) const;	// Checks whether a given index is inside of the domain.

private:
	int nX_ = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY_ = 0; // Amount of fields in the y-direction, not counting ghost cells.

	// todo: remove?
	inline int GetFlattenedIndex(const int xIdx, const int yIdx) const; // Helper function for getting the flattened index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].
	inline int GetFlattenedIndex(const CellIndex& cellIndex) const;
};