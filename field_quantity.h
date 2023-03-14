#include <vector>

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid)
struct FieldQuantity
{
	FieldQuantity()
	{};

	FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell=false, const double initialValue = 0.);

	//todo: Investigate how easy it is to switch the x/y order in which field vectors are stored, because this is going to cause a lot of headaches.

	// The value of the field at a specific x/y place. value[x][y]. the first column is the x coordinate, the second one the y coordinate, because of compatibility with legacy compatibility.
	std::vector<std::vector<double>> field;

	// Todo: implement the ghost cells here!

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value);
};