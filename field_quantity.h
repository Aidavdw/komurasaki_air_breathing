#include <vector>
#include <map>

enum EFieldQuantityBuffer
{
	MAIN, // The field that physically represents the value
	RUNGEKUTTA, // to be used for runge kutta iteration. A temporary buffer.
};

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid)
struct FieldQuantity
{
	FieldQuantity()
	{};

	FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const bool bCell=false, const double initialValue = 0.);

	//todo: Investigate how easy it is to switch the x/y order in which field vectors are stored, because this is going to cause a lot of headaches.

	// Todo: implement the ghost cells here!
	// The value of the field at a specific x/y place. value[x][y]. the first column is the x coordinate, the second one the y coordinate, because of compatibility with legacy compatibility.
	std::map<EFieldQuantityBuffer, std::vector<std::vector<double>>&> bufferMap;
	std::vector<std::vector<double>> field;
	std::vector<std::vector<double>> rungeKuttaBuffer;

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo=EFieldQuantityBuffer::MAIN);

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

private:
	// Sets a specific 2d vector to a value. Useful for internal use to call on the buffers.
	static void SetToValueInternal(std::vector<std::vector<double>>& valueField, const double value);
	static void Resize2DVector(std::vector<std::vector<double>>& valueField, const unsigned int sizeX, const unsigned int sizeY, const bool bCell, const double initialValue);
};