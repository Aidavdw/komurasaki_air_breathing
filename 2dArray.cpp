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
			// TODO: check if this actually calls the this() operator
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

void TwoDimensionalArray::ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to)
{
	// Needs to be a static function, because it's using the private data variable.
	for (size_t i = 0; i < (from.nX*from.nY); i++)
		to.data[i] = from.data[i];
}
