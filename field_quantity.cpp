#include "field_quantity.h"
#include <stdexcept>

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initialValue)
{
	// if it's set to be a cell, apparently according to florian it needs an additional element.
	Resize2DVector(field, sizeX, sizeY, bCell, initialValue);
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

}

void FieldQuantity::SetToValueInternal(std::vector<std::vector<double>>& Vec2D, const double value)
{
	for (int i = 0; i < Vec2D.size(); i++)
		for (int j = 0; j < Vec2D[0].size(); j++)
			Vec2D[i][j] = 0.0;
}

void FieldQuantity::Resize2DVector(std::vector<std::vector<double>>& vec2D, const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initialValue)
{
	// if it's set to be a cell, apparently according to florian it needs an additional element.
	std::vector<double> single_row(sizeX + bCell, initialValue);
	vec2D = std::vector<std::vector<double>>(sizeY + bCell, single_row);
}
