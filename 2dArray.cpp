#include "2dArray.h"
#include <stdexcept>

TwoDimensionalArray::TwoDimensionalArray(const unsigned int sizeX, const unsigned int sizeY, const double initialValue = 0, const double nGhostCells = 0) :
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

void TwoDimensionalArray::ElementWiseCopy(TwoDimensionalArray& from, TwoDimensionalArray& to)
{
	// Needs to be a static function, because it's using the private data variable.
	for (size_t i = 0; i < (from.nX*from.nY); i++)
		to.data[i] = from.data[i];
}
