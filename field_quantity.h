#include <vector>
#include <map>
#include "domain_enums.h"
#include "2dArray.h"

struct Domain;

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
	Domain* domain;

	int nGhostCells;

	TwoDimensionalArray main;				// The actual value of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray rungeKuttaBuffer;	// The Runge Kutta buffer of this field quantity. Access using At(), don't manually index!
	TwoDimensionalArray TBuffer;			// The T buffer of this field quantity. Access using At(), don't manually index!

	// Sets all the values in the field to this value.
	void SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo=EFieldQuantityBuffer::MAIN);

	void CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to);



	// operator overloaded accessor. Note that this is not the fastest way to set, so if possible do that directly on the 2d arrays level.
	inline double& operator () (int xIdx, int yIdx, EFieldQuantityBuffer buffer)
	{
		auto& buf = bufferMap.at(buffer);
		return buf[At(xIdx, yIdx)];
	}

	inline int At(const int xIdx, const int yIdx); // Helper function for getting the flattened index in the internal arrays for a certain index. Use like buffer_name[At(x,y)].

	inline int AtGhostCell(const EBoundaryLocation location, const int ghostX, const int ghostY); // Helper function for getting the flattened index of a ghost cell in the internal arrays for a certain index. Use like buffer_name[At(x,y)].

	double GetGradientInDirectionAndPosition(int posIdx[2], const double directionAngle);

private:
	int nX = 0; // Amount of fields in the x-direction, not counting ghost cells.
	int nY = 0; // Amount of fields in the y-direction, not counting ghost cells.
};