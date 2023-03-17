#include <vector>

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d arary.
struct TwoDimensionalArray
{
	TwoDimensionalArray()
	{};

	TwoDimensionalArray(const unsigned int sizeX, const unsigned int sizeY, const double initialValue = 0, const double nGhostCells = 0);

	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value);

	static void ElementWiseCopy(const TwoDimensionalArray& from, TwoDimensionalArray& to);

	void Resize(const unsigned int sizeX, const unsigned int sizeY, const double initialValue=0);

	bool IsEmpty() const;

	// operator overloaded accessor
	inline double& operator () (int xIdx, int yIdx)
	{
		// todo: check for the overhead of calling this with two numbers. Might be that the compiler doesnt pick up on it, and that [][] is actually faster.
		return data[(xIdx) + ((yIdx) * nX)];
	}

	// operator overloaded accessor
	inline double& operator [] (int flattenedIndex)
	{
		return data[flattenedIndex];
	}

private:
	std::vector<double> data;				// The actual value of this field quantity. Access using (), don't manually index!
};

