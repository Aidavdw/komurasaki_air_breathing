#include <vector>

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid)
struct FieldQuantity
{
	FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initalValue = 0.f);

	std::vector<std::vector<double>> field;
	

};