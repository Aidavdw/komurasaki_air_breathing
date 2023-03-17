#include "field_quantity.h"
#include <stdexcept>

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const double initialValue, const int nGhostCells = 2) : 
	nGhostCells(nGhostCells),
	nX(sizeX),
	nY(sizeY)
{
	main = TwoDimensionalArray(sizeX + nGhostCells, sizeY + nGhostCells, initialValue);
	rungeKuttaBuffer = TwoDimensionalArray(sizeX + nGhostCells, sizeY + nGhostCells, initialValue);
	TBuffer = TwoDimensionalArray(sizeX + nGhostCells, sizeY + nGhostCells, initialValue);

	bufferMap.insert({ EFieldQuantityBuffer::MAIN, main });
	bufferMap.insert({ EFieldQuantityBuffer::RUNGEKUTTA, rungeKuttaBuffer });
	bufferMap.insert({ EFieldQuantityBuffer::T, TBuffer });
}

void FieldQuantity::SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo = EFieldQuantityBuffer::MAIN)
{
	auto& buffer = bufferMap.at(bufferToWriteTo);
	buffer.SetAllToValue(value);
}

void FieldQuantity::CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to)
{
	auto& fromBuffer = bufferMap.at(from);
	auto& toBuffer = bufferMap.at(to);

	TwoDimensionalArray::ElementWiseCopy(fromBuffer, toBuffer);
}

inline int FieldQuantity::At(const int xIdx, const int yIdx)
{
	return (xIdx + nGhostCells) + ((yIdx + nGhostCells)*nX);
}

inline int FieldQuantity::AtGhostCell(const EBoundaryLocation location, const int ghostX, const int ghostY)
{
	int xOffset;
	int yOffset;

	switch (location)
	{
	case EBoundaryLocation::LEFT:
		xOffset = 0;
		yOffset = nGhostCells;
		break;
	case EBoundaryLocation::RIGHT:
		xOffset = nX;
		yOffset = nGhostCells;
		break;
	case EBoundaryLocation::TOP:
		xOffset = nGhostCells;
		yOffset = nY;
		break;
	case EBoundaryLocation::BOTTOM:
		xOffset = nGhostCells;
		yOffset = 0;
	}
	
	return (ghostX + xOffset) + ((ghostY + yOffset) * nX);
}