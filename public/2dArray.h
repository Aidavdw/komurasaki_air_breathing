#pragma once
#include <vector>
#include "domain_enums.h"
#include "index2d.h"

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d array.
class TwoDimensionalArray
{
public:
	TwoDimensionalArray() = default;

	TwoDimensionalArray(const int sizeX, const int sizeY, const int nGhostCells, const double initialValue = 0);

	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.
	int nGhostCells = 0; // The amount of ghostCells in this two-dimensional array. If zero, considered as just a 2d array.

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value);

	// Copies all the information from one array to another. Requires the arrays to be the same size.
	static void ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to);

	void Resize(const int sizeX, const int sizeY, const double initialValue=0);

	// Returns a transposed copy of this array.
	TwoDimensionalArray Transpose() const;

	// Returns whether or not all the entries are (approximately) zero. Note that this is an expensive operation, and should only be used for debugging, as it iterates over all the cells.
	bool IsFilledWithZeroes() const;

	bool IsEmpty() const;
	bool IsSquare() const;
	bool IsLowerTriangular() const;
	bool IsUpperTriangular() const;
	bool IsDiagonallySymmetric() const;
	bool HasDiagonalGrainsOnly(const int kernelSize) const; // Checks if the array only has nonzero entries on the diagonal, or a distance of kernelSize away from the diagonal.



	/* Getters for index */
	//todo: replace by CellIndex
	std::pair<int, int> GetIndexFromIndexOnBoundary(const EFace boundary, const int indexOnBoundary) const;
	
	// todo: check for the overhead of calling this with two numbers. Might be that the compiler doesnt pick up on it, and that [][] is actually faster.
	/* operator () overloads for get/setter: with pair of ints, cellIndex, and a special one for accessing ghostCells. */
	inline double& operator () (const int xIdx, const int yIdx)							{ return data_.at((xIdx + nGhostCells) + ((yIdx + 2*nGhostCells) * nX)); }
	inline double& operator () (const CellIndex& cellIndex)								{ return operator()(cellIndex.x, cellIndex.y); }
	inline double& GetReferenceIncludingGhostCells(const int xIdx, const int yIdx)		{ return data_.at((xIdx) + ((yIdx) * nX)); }
	inline double& GetReferenceIncludingGhostCells(const CellIndex& cellIndex)			{ return GetReferenceIncludingGhostCells(cellIndex.x, cellIndex.y); }

	/* operator () overloads for const getter: with pair of ints, cellIndex, and a special one for accessing ghostCells. */
	inline double GetAt(const int xIdx, const int yIdx) const					{ return data_.at((xIdx + nGhostCells) + ((yIdx + 2*nGhostCells) * nX)); }
	inline double GetAt(const CellIndex& cellIndex) const						{ return GetAt(cellIndex.x, cellIndex.y); }
	inline double GetIncludingGhostCells(const int xIdx, const int yIdx) const	{ return data_.at((xIdx) + ((yIdx) * nX)); }
	inline double GetIncludingGhostCells(const CellIndex& cellIndex) const		{ return GetIncludingGhostCells(cellIndex.x, cellIndex.y); }

private:
	std::vector<double> data_;				// The actual value of this field quantity. Access using (), don't manually index!
};

