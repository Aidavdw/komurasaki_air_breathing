#pragma once
#include <string>
#include <vector>
#include "domain_enums.h"
#include "index2d.h"
#include "muscl.h"

/*Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d array. If initialised with GhostCells, these are also included, but they are indexed relative to the original dimensions. This means that the left ghost cells will be -i, and the ones on the right x_tot + i (and vice versa for y).
 * Has const accessors GetAt() and setters At(), taking either int or CellIndex as an argument, with input validation for debug builds.
 */
class TwoDimensionalArray
{
public:
	TwoDimensionalArray() = default;

	TwoDimensionalArray(const int sizeX, const int sizeY, const int nGhostCells, const double initialValue = 0);

	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.
	int nGhostCells = 0; // The amount of ghostCells in this two-dimensional array. If zero, considered as just a 2d matrix.
	
	void SetAllToValue(const double value); // Sets all the values in the field to this value. Useful for initialisation or resetting.
	static void ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to); // Copies all the information from one array to another. Requires the arrays to be the same size.
	void Resize(const int sizeX, const int sizeY, const int setGhostCellCount = 0, const double initialValue = 0); // Changes the size of the array. Really only used for initialisation.
	TwoDimensionalArray Transpose() const; // Returns a transposed copy of this array.

	double GetMUSCLInterpolationForFace(const CellIndex& cix, const EFace face, const EAxisDirection faceNormalDirection, const double MUSCLBias, const EFluxLimiterType fluxLimiterType) const;
	
	bool IsFilledWithZeroes() const; // Returns whether or not all the entries are (approximately) zero. Note that this is an expensive operation, and should only be used for debugging, as it iterates over all the cells.
	bool IsEmpty() const;
	bool IsSquare() const;
	bool IsLowerTriangular() const;
	bool IsUpperTriangular() const;
	bool IsDiagonallySymmetric() const;
	bool HasDiagonalGrainsOnly(const int kernelSize) const; // Checks if the array only has nonzero entries on the diagonal, or a distance of kernelSize away from the diagonal.
	bool IsDiagonalBlockMatrix(const int blockSize) const;


	/******* INLINE ACCESSORS: reading off the value/slicing  ************/

	inline bool operator == (const TwoDimensionalArray& other) const;

	
	inline double& operator () (const CellIndex& cellIndex)						{ return operator()(cellIndex.x, cellIndex.y); }
	inline double GetAt(const CellIndex& cellIndex) const						{ return GetAt(cellIndex.x, cellIndex.y); }

	// todo: rename this to GetReferenceToGhostCell and flip the default bool, as this is only used to set the ghost cells.
	inline double& GetReferenceIncludingGhostCells(const CellIndex& cellIndex, const bool bAllowNormalCells=true)
	{
#ifdef _DEBUG
		if (cellIndex.relativeToBoundary != BOTTOM)
			throw std::logic_error("cellIndexes are defined relative to the bottom-left corner. A cellIndex based on a different boundary is supplied");
#endif
		return GetReferenceIncludingGhostCells(cellIndex.x, cellIndex.y, bAllowNormalCells);
	}
	
	inline double GetIncludingGhostCells(const CellIndex& cellIndex, const bool bAllowNormalCells=false) const
	{
#ifdef _DEBUG
		if (cellIndex.relativeToBoundary != BOTTOM)
			throw std::logic_error("cellIndexes are defined relative to the bottom-left corner. A cellIndex based on a different boundary is supplied");
#endif
		return GetIncludingGhostCells(cellIndex.x, cellIndex.y, bAllowNormalCells);
	}
	
	inline int GetRowOffset(const int xIdx, const int yIdx) const		{ return (yIdx + nGhostCells) * (2*nGhostCells + nX);	} // Total amount of elements to skip because of the amount of rows the precede it; there are two entire rows of ghost cells before index 0, so (yIdx + nGhostCells). Then, all the rows (including the skipped ones) are (2*nGhostCells + nX) long, because they have ghost cells on both sides.
	inline int GetColumnOffset(const int xIdx, const int yIdx) const	{ return xIdx + nGhostCells;	} // Amount of elements to skip in this row due to the ghost cells presence. Linearly accessed, so only skip 1*nGhostCells.
	
	/* operator () overloads for get/setter: with pair of ints, cellIndex, and a special one for accessing ghostCells. */
	inline double& operator () (const int xIdx, const int yIdx)
	{

		const int rowOffset = GetRowOffset(xIdx, yIdx);
		const int colOffset = GetColumnOffset(xIdx, yIdx);
#ifdef _DEBUG
		if (xIdx < 0 || xIdx >= nX || yIdx < 0 || yIdx >= nY)
			throw std::runtime_error("Tried accessing 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is out of range (or a ghost cell).");
		if (isnan(data_.at(rowOffset + colOffset)))
			throw std::runtime_error("Value at [" + std::to_string(xIdx) + ", " + std::to_string(yIdx) + "] is NaN" );
#endif
		return data_.at(rowOffset + colOffset);
	}

	// Setter for GhostCells
	inline double& GetReferenceIncludingGhostCells(const int xIdx, const int yIdx, const bool bAllowNormalCells=false)
	{
#ifdef  _DEBUG
		// if it's not in a ghost cell, throw error
		if (xIdx < -nGhostCells || xIdx >= nX + 2*nGhostCells || yIdx < -nGhostCells || yIdx >= nY + 2*nGhostCells)
			throw std::runtime_error("Tried accessing Ghost Cell of 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is out of range.");
		if (!bAllowNormalCells)
		{
			if ((xIdx > 0 && xIdx < nX) && (yIdx > 0 && yIdx < nY))
				throw std::runtime_error("Tried accessing 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is inside the normal domain. Be sure to use operator() to access these.");
		}
#endif
		const int rowOffset = GetRowOffset(xIdx, yIdx);
		const int colOffset = GetColumnOffset(xIdx, yIdx);
		return data_.at(rowOffset + colOffset);
	}

	/* operator () overloads for const getter: with pair of ints, cellIndex, and a special one for accessing ghostCells. */
	inline double GetAt(const int xIdx, const int yIdx) const
	{
		const int rowOffset = GetRowOffset(xIdx, yIdx);
		const int colOffset = GetColumnOffset(xIdx, yIdx);
#ifdef _DEBUG
		if (xIdx < 0 || xIdx >= nX || yIdx < 0 || yIdx >= nY)
			throw std::runtime_error("Tried accessing 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is out of range (or a ghost cell).");
		if (isnan(data_.at(rowOffset + colOffset)))
			throw std::runtime_error("Value at [" + std::to_string(xIdx) + ", " + std::to_string(yIdx) + "] is NaN" );
#endif
		return data_.at(rowOffset + colOffset);
	}

	// Getter for GhostCells
	inline double GetIncludingGhostCells(const int xIdx, const int yIdx, const bool bAllowNormalCells=true) const
	{
#ifdef  _DEBUG
		// if it's not in a ghost cell, throw error
		if (xIdx < -nGhostCells || xIdx >= nX + 2*nGhostCells || yIdx < -nGhostCells || yIdx >= nY + 2*nGhostCells)
			throw std::runtime_error("Tried accessing Ghost Cell of 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is out of range.");
		if (!bAllowNormalCells)
		{
			if ((xIdx > 0 && xIdx < nX) && (yIdx > 0 && yIdx < nY))
				throw std::runtime_error("Tried accessing 2D array at index [" + std::to_string(xIdx) + "," + std::to_string(yIdx) + "], which is inside the normal domain. Be sure to use operator() to access these.");
		}
#endif
		const int rowOffset = GetRowOffset(xIdx, yIdx);
		const int colOffset = GetColumnOffset(xIdx, yIdx);
		return data_.at(rowOffset + colOffset);
	}
	
	double GetAtPythonProxy(const int xIdx, const int yIdx) const { return GetAt(xIdx,yIdx); } // pybind11 does not like functions that take a lot of aliases, so here is a unique proxy.

private:
	std::vector<double> data_;				// The actual value of this field quantity. Access using (), don't manually index!
};

