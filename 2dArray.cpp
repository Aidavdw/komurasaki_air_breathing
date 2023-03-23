#include "2dArray.h"
#include <stdexcept>

TwoDimensionalArray::TwoDimensionalArray(const unsigned int sizeX, const unsigned int sizeY, const double initialValue = 0) :
	nX(sizeX),
	nY(sizeY)
{
	Resize(sizeX, sizeY, initialValue);
}

void TwoDimensionalArray::SetAllToValue(const double value)
{
	for (int i = 0; i < data.size(); i++)
		data[i] = value;
}

void TwoDimensionalArray::Resize(const unsigned int sizeX, const unsigned int sizeY, const double initialValue)
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

bool TwoDimensionalArray::IsEmpty() const
{
	if (nX == 0 || nY == 0)
		return true;
	return data.empty();
}

std::pair<int, int> TwoDimensionalArray::GetIndexFromIndexOnBoundary(const EBoundaryLocation boundary, const int indexOnBoundary) const
{
	switch (boundary)
	{
	case EBoundaryLocation::LEFT:
		return std::make_pair(0, indexOnBoundary);
	case EBoundaryLocation::RIGHT:
		return std::make_pair(nX - 1, indexOnBoundary);
	case EBoundaryLocation::BOTTOM:
		return std::make_pair(indexOnBoundary, 0);
	case EBoundaryLocation::TOP:
		return std::make_pair(indexOnBoundary, nY - 1);
	default:
		throw std::logic_error("FlattenIndexOnBoundary is not implemented for this boundary location.");
	}
}

void TwoDimensionalArray::ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to)
{
	if (from.nX != to.nX || from.nY != to.nY)
		throw std::logic_error("Cannot do element-wise copying between two TwoDimensionalArrays that are not the same size.");
	// Needs to be a static function, because it's using the private data variable.
	for (size_t i = 0; i < (from.nX*from.nY); i++)
		to.data[i] = from.data[i];
}
