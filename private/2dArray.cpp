#include "2dArray.h"
#include <stdexcept>
#include <cmath>

#include "AuxFunctions.h"

TwoDimensionalArray::TwoDimensionalArray(const int sizeX, const int sizeY, const double initialValue = 0) :
	nX(sizeX),
	nY(sizeY)
{
	Resize(sizeX, sizeY, initialValue);
}

void TwoDimensionalArray::SetAllToValue(const double value)
{
	for (double& i : data)
		i = value;
}

void TwoDimensionalArray::Resize(const int sizeX, const int sizeY, const double initialValue)
{
	data = std::vector<double>((sizeX) * (sizeY), initialValue);
}

TwoDimensionalArray TwoDimensionalArray::Transpose() const
{
	auto T = TwoDimensionalArray(nY, nX);

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
	for (const auto& val : data)
		sum += val;

	return (sum < 0.001);
}

bool TwoDimensionalArray::IsEmpty() const
{
	if (nX == 0 || nY == 0)
		return true;
	return data.empty();
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

std::pair<int, int> TwoDimensionalArray::GetIndexFromIndexOnBoundary(const EFace boundary, const int indexOnBoundary) const
{
	switch (boundary)
	{
	case EFace::LEFT:
		return std::make_pair(0, indexOnBoundary);
	case EFace::RIGHT:
		return std::make_pair(nX - 1, indexOnBoundary);
	case EFace::BOTTOM:
		return std::make_pair(indexOnBoundary, 0);
	case EFace::TOP:
		return std::make_pair(indexOnBoundary, nY - 1);
	default:
		throw std::logic_error("FlattenIndexOnBoundary is not implemented for this boundary location.");
	}
}

void TwoDimensionalArray::ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to)
{
	#ifdef _DEBUG
	if (from.nX != to.nX || from.nY != to.nY)
		throw std::logic_error("Cannot do element-wise copying between two TwoDimensionalArrays that are not the same size.");

	#endif
	
	// Needs to be a static function, because it's using the private data variable.
	for (size_t i = 0; i < (from.nX*from.nY); i++)
		to.data[i] = from.data[i];
}
