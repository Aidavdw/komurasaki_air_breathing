#include "field_quantity.h"

#include <cassert>

#include "pos2d.h"
#include "index2d.h"
#include "domain.h"
#include <stdexcept>
#include <cmath>

#include "muscl.h"

FieldQuantity::FieldQuantity(Domain* domain, const int sizeX, const int sizeY, const double initialValue, const int nGhostCells) :
	domain(domain),
	nGhostCells(nGhostCells),
	nX_(sizeX),
	nY_(sizeY)
{
	currentTimeStep = TwoDimensionalArray(sizeX + 2 * nGhostCells, sizeY + 2 * nGhostCells, initialValue);
	rungeKuttaBuffer = TwoDimensionalArray(sizeX + 2 * nGhostCells, sizeY + 2 * nGhostCells, initialValue);
	flux = TwoDimensionalArray(sizeX + 2 * nGhostCells, sizeY + 2 * nGhostCells, initialValue);

	bufferMap.insert({ EFieldQuantityBuffer::CURRENT_TIME_STEP, currentTimeStep });
	bufferMap.insert({ EFieldQuantityBuffer::RUNGE_KUTTA, rungeKuttaBuffer });
	bufferMap.insert({ EFieldQuantityBuffer::NEXT_TIME_STEP, nextTimeStepBuffer });
	bufferMap.insert({ EFieldQuantityBuffer::FLUX, flux });
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
double FieldQuantity::GetInterpolatedValueAtPosition(const Position& atPosition, const EFieldQuantityBuffer bufferName) const
{
	Position distanceFromCellCenter;
	const CellIndex cellIndex = domain->InvertPositionToIndex(atPosition, distanceFromCellCenter);

	// Setting which cells to interpolate with
	const CellIndex horizontalInterpolateTarget = (distanceFromCellCenter.x < 0) ? cellIndex + CellIndex{1,0} : cellIndex + CellIndex{-1, 0} ;
	const CellIndex verticalInterpolateTarget = (distanceFromCellCenter.y < 0) ? cellIndex + CellIndex{0,1} : cellIndex + CellIndex{0, -1} ;
	
	std::pair<double, double> cellSize = domain->GetCellSizes(cellIndex);
	std::pair<double, double> xInterpolateSize = domain->GetCellSizes(horizontalInterpolateTarget);
	std::pair<double, double> yInterpolateSize = domain->GetCellSizes(verticalInterpolateTarget);

	const TwoDimensionalArray& buffer = bufferMap.at(bufferName);
	double deltaHorizontal = buffer.GetAt(horizontalInterpolateTarget) + buffer.GetAt(cellIndex) / (0.5*(cellSize.first + xInterpolateSize.first)) - buffer.GetAt(cellIndex);
	double deltaVertical = buffer.GetAt(verticalInterpolateTarget) + buffer.GetAt(cellIndex) / (0.5*(cellSize.second + yInterpolateSize.second)) - buffer.GetAt(cellIndex);
	double interpolatedValue = buffer.GetAt(cellIndex) + deltaHorizontal + deltaVertical;
		
#ifdef _DEBUG
	// Check if it's between the value of cells that it's interpolating from as a sanity check.
	double highestVal = std::max({buffer.GetAt(horizontalInterpolateTarget), buffer.GetAt(verticalInterpolateTarget), buffer.GetAt(cellIndex)});
	assert(interpolatedValue < highestVal);
	double lowestVal = std::min({buffer.GetAt(horizontalInterpolateTarget), buffer.GetAt(verticalInterpolateTarget), buffer.GetAt(cellIndex)});
	assert(interpolatedValue > lowestVal);
#endif
	
	return interpolatedValue;
}

void FieldQuantity::PopulateMUSCLBuffers(const EFieldQuantityBuffer sourceBuffer, double MUSCLBias, const EFluxLimiterType fluxLimiterType)
{
	if (nGhostCells < 2)
		throw std::logic_error("Cannot populate MUSCL buffers, as the amount of ghost cells is smaller than two. Because of how MUSCL has been implemented, given that c is at the edge, the sampling at p2 will be done on the ghost cells! This means that there is a soft lower limit of 2 on the amount of ghost cells.");
	
	const TwoDimensionalArray& source = bufferMap.at(sourceBuffer);
	for (int xIdx = 0; xIdx < nX_; ++xIdx)
	{
		for (int yIdx = 0; yIdx < nY_; yIdx++)
		{
			
#ifdef _DEBUG
			// Check if the cells are not zero; that would probably be an error.
			if (IsCloseToZero(source.GetAt(xIdx-1, yIdx)))
				throw std::logic_error("MUSCL m cell is zero. Are you sure you're doing this right?");
			if (IsCloseToZero(source.GetAt(xIdx, yIdx)))
				throw std::logic_error("MUSCL centre cell is zero. Are you sure you're doing this right?");
			if (IsCloseToZero(source.GetAt(xIdx+1, yIdx)))
				throw std::logic_error("MUSCL p1 cell is zero. Are you sure you're doing this right?");
			if (IsCloseToZero(source.GetAt(xIdx, yIdx+2)))
				throw std::logic_error("MUSCL p2 cell is zero. Are you sure you're doing this right?");
#endif
			// Left/right is straightforward.
			// TODO: check if the left/right mixup here is correct.
			MUSCLBuffer[EFace::RIGHT](xIdx, yIdx) = MUSCLInterpolate(source.GetAt(xIdx-1, yIdx), source.GetAt(xIdx, yIdx), source.GetAt(xIdx+1, yIdx), source.GetAt(xIdx+2, yIdx), EMUSCLSide::LEFT, MUSCLBias, fluxLimiterType);
			MUSCLBuffer[EFace::LEFT](xIdx, yIdx) = MUSCLInterpolate(source.GetAt(xIdx-1, yIdx), source.GetAt(xIdx, yIdx), source.GetAt(xIdx+1, yIdx), source.GetAt(xIdx+2, yIdx), EMUSCLSide::RIGHT, MUSCLBias, fluxLimiterType);

			// MUSCL ON TOP FACE -> Top face flux (Left = Down and Right = Up)
			MUSCLBuffer[EFace::BOTTOM](xIdx, yIdx) = MUSCLInterpolate(source.GetAt(xIdx, yIdx-1), source.GetAt(xIdx, yIdx), source.GetAt(xIdx, yIdx+1), source.GetAt(xIdx, yIdx+2), EMUSCLSide::LEFT, MUSCLBias, fluxLimiterType);
			MUSCLBuffer[EFace::TOP](xIdx, yIdx) = MUSCLInterpolate(source.GetAt(xIdx, yIdx-1), source.GetAt(xIdx, yIdx), source.GetAt(xIdx, yIdx+1), source.GetAt(xIdx, yIdx+2), EMUSCLSide::RIGHT, MUSCLBias, fluxLimiterType);
		}
	}
	
}

double& FieldQuantity::operator()(const int xIdx, const int yIdx, const EFieldQuantityBuffer buffer)
{
#ifdef _DEBUG
	const CellIndex c = {xIdx, yIdx};
	if (!IsValidIndex(c, false))
		throw std::logic_error("Cannot access field quantity at index" + c.ToString());
#endif
		
	auto& buf = bufferMap.at(buffer);
	return buf[GetFlattenedIndex(xIdx, yIdx)];
}

inline int FieldQuantity::GetFlattenedIndex(const int xIdx, const int yIdx) const
{
	return (xIdx + nGhostCells) + ((yIdx + nGhostCells)*nX_);
}
int FieldQuantity::GetFlattenedIndex(const CellIndex &cellIndex) const
{
	return (cellIndex.x + nGhostCells) + ((cellIndex.y + nGhostCells)*nX_);
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
	double partialDerivativeX = (currentTimeStep.GetAt(rootPos.x, rootPos.y) - currentTimeStep.GetAt(dxPos.x, dxPos.y)) / dx;
	double partialDerivativeY = (currentTimeStep.GetAt(rootPos.x, rootPos.y) - currentTimeStep.GetAt(dyPos.x, dyPos.y)) / dy;

	// Turning the direction given into a unit vector, u = < cos(theta), sin(theta) >. Since this has unit length, decomposing the partial derivatives calculated above is really straightforward.
	return cos(directionAngle) * partialDerivativeX + sin(directionAngle) * partialDerivativeY;

}
double FieldQuantity::GetAverageValue(const bool bExpectUniformField) const
{
	#ifdef _DEBUG
	if (bExpectUniformField)
	{
		// Even if it's a uniform field, just do a littttttle check to make sure it actually is
		double checkval1 = currentTimeStep.GetAt(currentTimeStep.nX/4, currentTimeStep.nY/2);
		double checkval2 = currentTimeStep.GetAt(3*currentTimeStep.nX/4, currentTimeStep.nY/2);
		if (checkval1 - checkval2 > 0.005)
			throw std::logic_error("Expected uniform field, but it was not!");
	}
	#endif

	// sort of arbitrary where it reads from.
	if (bExpectUniformField)
		return currentTimeStep.GetAt(currentTimeStep.nX/4, currentTimeStep.nY/2);

	double totalSum = 0;
	for (int xIdx = 0; xIdx < nX_; xIdx++)
	{
		for (int yIdx = 0; yIdx < nY_; yIdx++)
		{
			totalSum += currentTimeStep.GetAt(xIdx,yIdx);
		}
	}
	
	return totalSum / (nX_ * nY_);
}

bool FieldQuantity::IsValidIndex(const CellIndex& cellIndex, const bool bAllowGhostCells) const
{
	if (cellIndex.x < -nGhostCells*bAllowGhostCells)
		return false;
	if (cellIndex.x > nX_ + nGhostCells*bAllowGhostCells)
		return false;
	if (cellIndex.y < -nGhostCells*bAllowGhostCells)
		return false;
	if (cellIndex.y > nY_ + nGhostCells*bAllowGhostCells)
		return false;
	
	return true;
}

