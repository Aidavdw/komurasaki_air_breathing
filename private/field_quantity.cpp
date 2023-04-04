#include "field_quantity.h"
#include "domain.h"
#include <stdexcept>
#include <cmath>

FieldQuantity::FieldQuantity(Domain* domain, const int sizeX, const int sizeY, const double initialValue, const int nGhostCells) :
	domain(domain),
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

void FieldQuantity::SetAllToValue(const double value, const EFieldQuantityBuffer bufferToWriteTo)
{
	auto& buffer = bufferMap.at(bufferToWriteTo);
	buffer.SetAllToValue(value);
}

void FieldQuantity::CopyToBuffer(const EFieldQuantityBuffer from, const EFieldQuantityBuffer to)
{
	TwoDimensionalArray::ElementWiseCopy(bufferMap.at(from), bufferMap.at(to));
}

inline int FieldQuantity::At(const int xIdx, const int yIdx) const
{
	return (xIdx + nGhostCells) + ((yIdx + nGhostCells)*nX);
}

inline int FieldQuantity::AtGhostCell(const EBoundaryLocation location, const int ghostX, const int ghostY) const
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
		break;
	default:
		throw std::logic_error("Can't check a ghost cell from this side.");
	}
	
	return (ghostX + xOffset) + ((ghostY + yOffset) * nX);
}

double FieldQuantity::GetGradientInDirectionAndPosition(const CellIndex posIdx, const double directionAngle) const
{
	// The angle will be decomposed in an x- component and a y component. Using a forward difference scheme, the two will then be linearly interpolated to get a gradient in the direction angle (theta) direction.

	// Could do higher order finite difference method too?

	// todo: Right now, the entire cell is assumed to have the same value over its entire area, with a jump at the edge. The nearest cell center is considered only. Possible expansion would be to interpolate the values at these delta poses based on their neighbours.

	/**** First get the partial derivatives in the xand y directions. ***/
	CellIndex rootPos = posIdx;
	CellIndex dxPos = { posIdx.x - 1 , posIdx.y };
	CellIndex dyPos = { posIdx.x, posIdx.y - 1};
	
	// can't get values outside of the domain, so clamp if this is at the edge. Backward difference.
	if (posIdx.x == 0)
	{
		rootPos.x += 1;
		dxPos.x = 0;
	}
	else if (posIdx.y == 0)
	{
		dyPos.y += 0;
		rootPos.y += 1;
	}

	// Calculating the partial derivatives in the x-y directions.
	double dx = (domain->meshSpacing[0].GetCellWidth(dxPos.x) + domain->meshSpacing[0].GetCellWidth(rootPos.x)) * 0.5;
	double dy = (domain->meshSpacing[1].GetCellWidth(dyPos.y) + domain->meshSpacing[1].GetCellWidth(rootPos.y)) * 0.5;
	double partialDerivativeX = (main.GetAt(rootPos.x, rootPos.y) - main.GetAt(dxPos.x, dxPos.y)) / dx;
	double partialDerivativeY = (main.GetAt(rootPos.x, rootPos.y) - main.GetAt(dyPos.x, dyPos.y)) / dy;

	// Turning the direction given into a unit vector, u = < cos(theta), sin(theta) >. Since this has unit length, decomposing the partial derivatives calculated above is really straightforward.
	return cos(directionAngle) * partialDerivativeX + sin(directionAngle) * partialDerivativeY;

}
