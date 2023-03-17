#include <vector>
#include <map>
#include "domain_enums.h"
#include "2dArray.h"

enum EFieldQuantityBuffer
{
	MAIN,			// The field that physically represents the value
	RUNGEKUTTA,		// to be used for runge kutta iteration. A temporary buffer.
	T,				// not yet sure what this means but oh well!
};

// Thin wrapper around a 2d array with doubles, representing a value over an entire domain (in a grid). It is implemented as a flattened 2d arary.
struct FieldQuantity
{
	FieldQuantity() :
		nGhostCells(0)
	{};

	FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const double initialValue = 0., const int nGhostCells=2);

	std::map<EFieldQuantityBuffer, TwoDimensionalArray> bufferMap;

	int nGhostCells;

	TwoDimensionalArray main;				// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;	// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray TBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo=EFieldQuantityBuffer::MAIN);

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);

	inline int At(const int xIdx, const int yIdx); // Helper function for getting the index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].
	inline int AtGhostCell(const EBoundaryLocation location, const int ghostX, const int ghostY); // Helper function for getting a ghost cell in the internal arrays for a certain index.Use like buffer_name[At(boundaryLocation, x, y)].

private:
	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.
};