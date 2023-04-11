#pragma once
#include <vector>
#include "domain_enums.h"
#include "index2d.h"

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d array.
struct TwoDimensionalArray
{
	TwoDimensionalArray() = default;

	TwoDimensionalArray(const int sizeX, const int sizeY, const double initialValue = 0);

	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value);

	// Copies all the information from one array to another. Requires the arrays to be the same size.
	static void ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to);

	void Resize(const int sizeX, const int sizeY, const double initialValue=0);

	// Returns a transposed copy of this array.
	TwoDimensionalArray Transpose() const;

	bool IsEmpty() const;

	bool IsSquare() const;
	bool IsLowerTriangular() const;
	bool IsUpperTriangular() const;
	bool IsSymmetric() const;



	/* Getters for index */
	std::pair<int, int> GetIndexFromIndexOnBoundary(const EBoundaryLocation boundary, const int indexOnBoundary) const;

	/* Getters for actual values */

	// operator overloaded accessor. This can be used for setting values.
	inline double& operator () (int xIdx, int yIdx)
	{
		// todo: check for the overhead of calling this with two numbers. Might be that the compiler doesnt pick up on it, and that [][] is actually faster.
		return data[(xIdx) + ((yIdx) * nX)];
	}

	// operator overloaded accessor. This can be used for setting values.
	inline double& operator () (const CellIndex& cellIndex)
	{
		// todo: check for the overhead of calling this with two numbers. Might be that the compiler doesnt pick up on it, and that [][] is actually faster.
		return data[(cellIndex.x)+((cellIndex.y)*nX)];
	}

	// const operator overloaded getter. Note that does does not allow setting.
	inline double GetAt(int xIdx, int yIdx) const
	{
		return data[(xIdx)+((yIdx)*nX)];
	}

	// operator overloaded accessor
	inline double& operator [] (int flattenedIndex)
	{
		return data[flattenedIndex];
	}

private:
	std::vector<double> data;				// The actual value of this field quantity. Access using (), don't manually index!
};

