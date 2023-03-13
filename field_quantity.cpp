#include "field_quantity.h"

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initialValue)
{
	// if it's set to be a cell, apparently according to florian it needs an additional element.
	std::vector<double> single_row(sizeX + bCell, 0);
	field = std::vector<std::vector<double>>(sizeY + bCell, single_row);

	SetAllToValue(initialValue);
}

void FieldQuantity::SetAllToValue(const double value)
{
	for (int i = 0; i < field.size(); i++)
		for (int j = 0; j < field[0].size(); j++)
			field[i][j] = 0.0;
}
