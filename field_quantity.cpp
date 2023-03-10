#include "field_quantity.h"

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initalValue)
{
	// if it's set to be a cell, apparently according to florian it needs an additional element.
	std::vector<double> single_row(sizeX + bCell, 0);
	field = std::vector<std::vector<double>>(sizeY + bCell, single_row);

	for (int i = 0; i < field.size(); i++)
		for (int j = 0; j < single_row.size(); j++)
			field[i][j] = 0.0;
}
