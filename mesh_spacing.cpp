#include "mesh_spacing.h"
#include "field_quantity.h"
#include "domain.h"
#include <stdexcept>

MeshSpacing::MeshSpacing(const EMeshSpacingType meshSpacingType, const double elementLength, const int nElements, const double resolution_left, const double resolution_center, const double resolution_right) :
	spacingType(meshSpacingType),
	length(elementLength),
	amountOfElements(nElements)
{

	// As not all spacing types require all the spacing parameters to be filled in, it needs to be handled separately for the mesh spacing type. Still fill in all the other parameters so that they can be used for debugging and/or inspection.
	FitSpacingToParameters();
}

double MeshSpacing::GetCellWith(const int i)
{
	switch(spacingType)
	{
	case EMeshSpacingType::CONSTANT:
		return left;
	case EMeshSpacingType::LINEAR:
		return left + (right - left) * (double(i) / (amountOfElements - 1));
	default:
		throw std::logic_error("Getting cell width for this spacing type has not yet been implemented.");
	}
}

void MeshSpacing::FitSpacingToParameters()
{
	double origLeft = left;
	double origRight = right;
	double origCenter = center;
	int origAmountOfElements = amountOfElements;

	int degreesOfFreedom = IsCloseToZero(left) + IsCloseToZero(right) + IsCloseToZero(center) + (amountOfElements == 0);

	int requiredParameterCount;
	switch (spacingType)
	{
	case EMeshSpacingType::CONSTANT:
		requiredParameterCount = 1;
		break;
	case EMeshSpacingType::LINEAR:
		requiredParameterCount = 2;
	default:
		throw std::logic_error("Fitting this type of mesh spacing to parameters has not yet been implemented.");
	}

	if (degreesOfFreedom > 4 - requiredParameterCount)
		throw std::logic_error("The spacing type is underconstrained! More parameters need to be set.");
	if (degreesOfFreedom < 4 - requiredParameterCount)
		throw std::logic_error("The spacing type is overconstrained! More parameters need to be set.");

	// Do the actual gradient descent here
}

bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() )
{
	return std::abs(x) < tolerance;
}


