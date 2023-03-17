#include <vector>

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d arary.
struct TwoDimensionalArray
{
	TwoDimensionalArray()
	{};

	TwoDimensionalArray(const unsigned int sizeX, const unsigned int sizeY, const double initialValue = 0);



	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value);

	static void ElementWiseCopy(TwoDimensionalArray& from, TwoDimensionalArray& to);

	// operator overloaded accessor
	inline double& operator () (int xIdx, int yIdx)
	{
		return data[(xIdx) + ((yIdx)*nX)];
	}

private:
	std::vector<double> data;				// The actual value of this field quantity. Access using (), don't manually index!

	static void Resize(std::vector<double>& valueField, const unsigned int sizeX, const unsigned int sizeY, const double initialValue);

};

