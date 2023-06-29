#include "2dArray.h"

#include <assert.h>
#include <stdexcept>
#include <cmath>

#include "AuxFunctions.h"

TwoDimensionalArray::TwoDimensionalArray(const int sizeX, const int sizeY, const int nGhostCells, const double initialValue) :
	nX(sizeX),
	nY(sizeY),
	nGhostCells(nGhostCells)
{
	Resize(sizeX + 2*nGhostCells, sizeY + 2*nGhostCells, initialValue);
}

void TwoDimensionalArray::SetAllToValue(const double value)
{
	for (double& i : data_)
		i = value;
}

void TwoDimensionalArray::Resize(const int sizeX, const int sizeY, const double initialValue)
{
	if (sizeX == nX && sizeY == nY)
	{
		SetAllToValue(initialValue);
		return;
	}
	data_ = std::vector<double>((sizeX) * (sizeY), initialValue);
	nX = sizeX;
	nY = sizeY;
}

TwoDimensionalArray TwoDimensionalArray::Transpose() const
{
	auto T = TwoDimensionalArray(nY, nX, nGhostCells);

	for (int xIdx = 0; xIdx < nX; xIdx++)
	{
		for (int yIdx = 0; yIdx < nY; yIdx++)
		{
			T(yIdx, xIdx) = this->GetAt(xIdx, yIdx);
		}
	}

	return T;
}

bool TwoDimensionalArray::IsFilledWithZeroes() const
{
	double sum = 0;
	for (const auto& val : data_)
		sum += val;

	return (sum < 0.001);
}

bool TwoDimensionalArray::IsEmpty() const
{
	if (nX == 0 || nY == 0)
		return true;
	return data_.empty();
}
bool TwoDimensionalArray::IsSquare() const
{
	return (nX == nY);
}

bool TwoDimensionalArray::IsLowerTriangular() const
{
	if (!IsSquare())
		throw std::logic_error("Array is not square.");
	
	for (int yIndex = 0; yIndex < nY; yIndex++)
	{
		for (int xIndex = yIndex; xIndex < nX; xIndex++)
		{
			if (IsCloseToZero(GetAt(xIndex, yIndex)))
				return false;
		}
	}
	return true;
}

bool TwoDimensionalArray::IsUpperTriangular() const
{
	if (!IsSquare())
		throw std::logic_error("Array is not square.");
	
	for (int yIndex = 0; yIndex < nX; yIndex++)
	{
		for (int xIndex = 0; xIndex < yIndex; xIndex++)
		{
			if (IsCloseToZero(GetAt(xIndex, yIndex)))
				return false;
		}
	}
	return true;
}

bool TwoDimensionalArray::IsDiagonallySymmetric() const
{
	// An example of a diagonally symmetric matrix:
	/*
	 * X 1 2 3
	 * 1 X 4 5
	 * 2 4 X 6
	 * 3 5 6 X
	 */
	
	if (!IsSquare())
		throw std::logic_error("Array is not square.");
	
	for (int yIndex = 0; yIndex < nX; yIndex++)
	{
		for (int xIndex = 0; xIndex < yIndex; xIndex++)
		{
			if (IsCloseToZero(GetAt(xIndex, yIndex) - GetAt(yIndex, xIndex)))
				return false;
		}
	}
	return true;
	
}

bool TwoDimensionalArray::HasDiagonalGrainsOnly(const int kernelSize) const
{	
	// An example of a diagonally matrix with kernel size 2
	/*
	 * 1 1 0 0
	 * 1 1 1 0
	 * 0 1 1 1
	 * 0 0 1 1
	 */

	// An example of a diagonally matrix with kernel size 4
	/*
	 * 1 1 1 1 0 0 0 0 0 0
	 * 1 1 1 1 0 0 0 0 0 0
	 * 1 1 2 2 2 2 0 0 0 0
	 * 1 1 2 2 2 2 0 0 0 0
	 * 0 0 2 2 3 3 3 3 0 0
	 * 0 0 2 2 3 3 3 3 0 0
	 * 0 0 0 0 3 3 4 4 4 4
	 * 0 0 0 0 3 3 4 4 4 4
	 * 0 0 0 0 0 0 4 4 4 4
	 * 0 0 0 0 0 0 4 4 4 4
	 */

	// Can only work if the dimension is a multiple of the kernel size
	if (!IsSquare())
		throw std::logic_error("Array is not square.");
	if (2*nX % kernelSize != 0)
		throw std::logic_error("Cannot check if the matrix consists of grains only, as the matrix size is not a multiple of (half) the kernel size.");
	if (kernelSize % 2 != 0)
		throw std::logic_error("Diagonal grains check only works for even numbered grain sizes!");


	const int amountOfGrains = (2*nX/kernelSize);
	for (int yIndex = 0; yIndex < nX; yIndex++)
	{
		// Check if we're inside the kernel
		const int offsetNumber = static_cast<int>(kernelSize/2) * static_cast<int>(std::floor(2*yIndex/kernelSize)); // little bit of a magic function, please check the drawings above!
			
		int grainStartIndex = kernelSize + offsetNumber;
		int grainEndIndex = yIndex - kernelSize / 2 + offsetNumber;
		
		for (int xIndex = 0; xIndex < yIndex; xIndex++)
		{
			if (xIndex < grainStartIndex || xIndex > grainEndIndex)
				if (!IsCloseToZero(GetAt(xIndex,yIndex)))
					return false;
		}
	}
	return true;
}

bool TwoDimensionalArray::IsDiagonalBlockMatrix(const int blockSize) const
{
	/* looks like (blockSize=2)
	 * 1 1 0 0 0 0
	 * 1 1 0 0 0 0
	 * 0 0 1 1 0 0
	 * 0 0 1 1 0 0
	 * 0 0 0 0 1 1
	 * 0 0 0 0 1 1
	 */
	if (!IsSquare())
		throw std::logic_error("Array is not square.");
	if (nX % blockSize != 0)
		return false;
	

	// Create a copy of just the ones in the block places, then check if they're equal.
	auto blockOnlyCopy = TwoDimensionalArray(nX, nY, nGhostCells);
	int nBlocks = nX/blockSize;
	for (int blockIdx = 0; blockIdx < nBlocks; blockIdx++)
	{
		for (int i = 0; i < blockSize; i++)
		{
			for (int j = 0; j < blockSize; j++)
			{
				blockOnlyCopy(blockSize*blockIdx + i, blockSize*blockIdx + j) = GetAt(blockSize*blockIdx + i, blockSize*blockIdx + j);
			}
		}
	}
	return (*this == blockOnlyCopy);
}

bool TwoDimensionalArray::operator==(const TwoDimensionalArray& other) const
{
	if (nX != other.nX || nY != other.nY || nGhostCells != other.nGhostCells)
		return false;

	for (int i = 0; i < nX; i++)
	{
		for (int j = 0; j < nY; j++)
		{
			if (!IsCloseToZero(GetAt(i, j) - other.GetAt(i,j)))
				return false;
		}
	}
	return true;
}

void TwoDimensionalArray::ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to)
{
	if (from.nX != to.nX || from.nY != to.nY)
		throw std::logic_error("Cannot do element-wise copying between two TwoDimensionalArrays that are not the same size.");
	assert(from.nGhostCells == to.nGhostCells);
	
	// Needs to be a static function, because it's using the private data variable.
	for (size_t i = 0; i < from.data_.size(); i++) // Also does the ghost cells. Can be excluded, but more cluttered code.
		to.data_[i] = from.data_[i];
}
