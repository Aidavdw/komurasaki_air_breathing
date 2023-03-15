#include "field_quantity.h"
#include <stdexcept>

FieldQuantity::FieldQuantity(const unsigned int sizeX, const unsigned int sizeY, const double initialValue, const int nGhostCells = 2) : 
	nGhostCells(nGhostCells),
	nX(sizeX),
	nY(sizeY)
{
	Resize2DVector(main, sizeX, sizeY, initialValue, nGhostCells);

	bufferMap.insert({ EFieldQuantityBuffer::MAIN, main });
	bufferMap.insert({ EFieldQuantityBuffer::RUNGEKUTTA, rungeKuttaBuffer });
	bufferMap.insert({ EFieldQuantityBuffer::T, TBuffer });
}

void FieldQuantity::SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo = EFieldQuantityBuffer::MAIN)
{
	auto& buffer = bufferMap.at(bufferToWriteTo);
	SetToValueInternal(buffer, value);
}

void FieldQuantity::CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to)
{
	auto& fromBuffer = bufferMap.at(from);
	auto& toBuffer = bufferMap.at(to);

	for (size_t i = 0; i < fromBuffer.size(); i++)
		toBuffer[i] = fromBuffer[i];
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

void FieldQuantity::SetToValueInternal(std::vector<double>& Vec2D, const double value)
{
	// The ghost cells can be skipped, saves a tiny bit of time? For now not implemented.
	for (int i = 0; i < Vec2D.size(); i++)
		Vec2D[i] = 0.0;
}

void FieldQuantity::Resize2DVector(std::vector<double>& vec2D, const unsigned int sizeX, const unsigned int sizeY, const double initialValue, const int nGhostCells)
{
	vec2D = std::vector<double>( (sizeX+ 2*nGhostCells) * (sizeY + 2*nGhostCells), initialValue);
}
