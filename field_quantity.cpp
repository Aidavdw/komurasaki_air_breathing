#include "field_quantity.h"
#include "domain.h"
#include <stdexcept>
#include <cmath>

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

double FieldQuantity::GetGradientInDirectionAndPosition(int posIdx[2], const double directionAngle)
{
	// The angle will be decomposed in an x- component and a y component. Using a forward difference scheme, the two will then be linearly interpolated to get a gradient in the direction angle (theta) direction.

	// Could do higher order finite difference method too?
	

	/**** First get the partial derivatives in the xand y directions. ***/
	int dxPos[2] = { posIdx[0] - 1 , posIdx[1] };
	int dyPos[2] = { posIdx[0], posIdx[1] - 1};

	// can't get values outside of the domain, so clamp if this is at the edge.
	if (posIdx[0] == 0)
	{
		posIdx[0] += 1;
		dxPos[0] = 0;
	}
	else if (posIdx[1] == 0)
	{
		dyPos[0] += 0;
		posIdx[0] += 1;
	}

	double dx = (domain->meshSpacing[0].GetCellWidth(dxPos[0]) + domain->meshSpacing[0].GetCellWidth(posIdx[0])) * 0.5;
	double dy = (domain->meshSpacing[1].GetCellWidth(dyPos[1]) + domain->meshSpacing[1].GetCellWidth(posIdx[1])) * 0.5;

	
	double partialDerivativeX = (main.GetAt(posIdx[0], posIdx[1]) - main.GetAt(dxPos[0], dxPos[1])) / dx;
	double partialDerivativeY = (main.GetAt(posIdx[0], posIdx[1]) - main.GetAt(dyPos[0], dyPos[1])) / dy;

	/***** Project the given angle onto these two partial derivatives ****/
	return cos(directionAngle) * partialDerivativeX + sin(directionAngle) * partialDerivativeY;

}
