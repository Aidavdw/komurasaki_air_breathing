#pragma once
#include <map>
#include "domain_enums.h"
#include "2dArray.h"
#include "index2d.h"

struct Position;
struct Domain;

enum EFieldQuantityBuffer
{
	MAIN,			// The field that physically represents the value
	RUNGEKUTTA,		// to be used for runge kutta iteration. A temporary buffer.
	NEXTITER,				// not yet sure what this means but oh well!
};

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d arary.
struct FieldQuantity
{
	FieldQuantity() :
		domain(nullptr),
		nGhostCells(0)
	{}

	FieldQuantity(Domain* domain, const int sizeX, const int sizeY, const double initialValue = 0., const int nGhostCells=2);

	std::map<EFieldQuantityBuffer, TwoDimensionalArray> bufferMap;
	Domain* domain;

	int nGhostCells;

	TwoDimensionalArray main;				// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;	// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray TBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo=EFieldQuantityBuffer::MAIN);

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	// Instead of taking considering that the entire volume has one value (lumped parameter estimation), consider that the cell center has the value, and anything outside of that will be linearly interpolated between those cells.
	double GetInterpolatedValueAtPosition(const Position atPosition) const;



	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const int xIdx, const int yIdx, const EFieldQuantityBuffer buffer)
	{
		auto& buf = bufferMap.at(buffer);
		return buf[At(xIdx, yIdx)];
	}

	// Shorthand overloaded accessor for main buffer.
	inline double& operator () (const int xIdx, const int yIdx)
	{
		return main[At(xIdx, yIdx)];
	}

	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (const CellIndex& cellIndex, const EFieldQuantityBuffer buffer)
	{
		auto& buf = bufferMap.at(buffer);
		return buf[At(cellIndex.x, cellIndex.y)];
	}

	// Shorthand overloaded accessor for main buffer.
	inline double& operator () (const CellIndex& cellIndex)
	{
		return main[At(cellIndex.x, cellIndex.y)];
	}

	inline int At(const int xIdx, const int yIdx) const; // Helper function for getting the flattened index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].
	inline int At(const CellIndex& cellIndex) const;

	// Helper function for getting the flattened index of a ghost cell in the internal arrays for a certain index. ghostX and ghostY are in boundary-relative coordinate frame.Use like buffer_name[AtGhostCell(x,y)].
	inline Position AtGhostCell(const EBoundaryLocation location, Position& posInBoundaryReferenceFrame) const; 
	double GetGradientInDirectionAndPosition(const CellIndex posIdx, const double directionAngle) const;

	double GetAverageValue(const bool bExpectUniformField) const;

private:
	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.
};