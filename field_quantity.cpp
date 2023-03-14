#include "field_quantity.h"
#include <stdexcept>

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const double initialValue)
{
	Resize2DVector(field, sizeX, sizeY, initialValue);

	columnSize = sizeY;
}

void FieldQuantity::SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo = EFieldQuantityBuffer::MAIN)
{
	auto& buffer = bufferMap.at(bufferToWriteTo);
	SetToValueInternal(buffer, value);
}

void FieldQuantity::CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to)
{
	auto& fromBuffer = bufferMap.at(from);
	auto& toBuffer = bufferMap.at(to);

	for (size_t i = 0; i < fromBuffer.size(); i++)
		toBuffer[i] = fromBuffer[i];
}

inline int FieldQuantity::At(const int x, const int y)
{
	return x + (y*columnSize);
}

void FieldQuantity::SetToValueInternal(std::vector<double>& Vec2D, const double value)
{
	for (int i = 0; i < Vec2D.size(); i++)
		Vec2D[i] = 0.0;
}

void FieldQuantity::Resize2DVector(std::vector<double>& vec2D, const unsigned int sizeX, const unsigned int sizeY, const double initialValue)
{
	vec2D = std::vector<double>(sizeX * sizeY, initialValue);
}
