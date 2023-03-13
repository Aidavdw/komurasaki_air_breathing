#include "mesh_spacing.h"
#include <stdexcept>

MeshSpacing::MeshSpacing(const EMeshSpacingType meshSpacingType, const double elementLength, const int nElements, const double resolution_left, const double resolution_center, const double resolution_right) :
	spacingType(meshSpacingType),
	length(elementLength),
	amountOfElements(nElements)
{

	// As not all spacing types require all the spacing parameters to be filled in, it needs to be handled separately for the mesh spacing type. Still fill in all the other parameters so that they can be used for debugging and/or inspection.
	switch (meshSpacingType)
	{
	case EMeshSpacingType::CONSTANT:
		// Ignore the values that are used for resolutions entirely.
		left = length / amountOfElements;
		center = length / amountOfElements;
		right = length / amountOfElements;
		break;

	case EMeshSpacingType::LINEAR:
		// Can only have 1 be non-zero, otherwise the problem is overconstrained.
		int bNonZeroMembers = IsCloseToZero(resolution_left) + IsCloseToZero(resolution_left);
		if (!IsCloseToZero(resolution_center) && !bNonZeroMembers)
			throw std::invalid_argument("A linear mesh spacing cannot take only a center argument, as this will result in an underconstrained problem. Please supply only a left- or right argument.");
		if (bNonZeroMembers == 0 || bNonZeroMembers > 1)
			throw std::invalid_argument("A linear mesh spacing must take 1 non-zero spacing parameter.");

		// Depending on which mesh spacing parameter was given, fill in the other ones considering a linear distribution. Note that in reality, only the left value is used to calculate the value.
		if (!IsCloseToZero(resolution_left))
		{
			left = resolution_left;
			//todo: check this calculation, i dont like this integration
			right = length / (amountOfElements * resolution_left) + resolution_left;
			center = left + (right - left) / 2;
		}
		else if (!IsCloseToZero(resolution_center))
		{
			throw std::invalid_argument("A linear mesh spacing cannot take only a center argument, as this will result in an underconstrained problem. Please supply only a left- or right argument.");
		}
		else if (!IsCloseToZero(resolution_right))
		{
			right = resolution_right;
			//todo: check this calculation, i dont like this integration
			left = resolution_right - length / (amountOfElements * resolution_right);
			center = left + (right - left) / 2;
		}
		break;

	case EMeshSpacingType::PARABOLIC:
		throw std::logic_error("parabolic mesh spacing is not yet implemented.");
	case EMeshSpacingType::EXPONENTIAL:
		throw std::logic_error("exponential mesh spacing is not yet implemented.");
	default:
		throw std::logic_error("Invalid Mesh spacing type.");
	}
}

double MeshSpacing::GetCellSize(const int index)
{
	switch (spacingType)
	{
	case EMeshSpacingType::CONSTANT:
		return length / amountOfElements;
	case EMeshSpacingType::LINEAR:
		double nondimensionalised_x = index - amountOfElements / double(amountOfElements);
		return nondimensionalised_x * (right - left);
	case EMeshSpacingType::PARABOLIC:
		throw std::logic_error("parabolic mesh spacing is not yet implemented.");
	case EMeshSpacingType::EXPONENTIAL:
		throw std::logic_error("exponential mesh spacing is not yet implemented.");
	}
	
	return 0;
}

double MeshSpacing::GetCellPosition(const int index)
{
	return 0.0;
}

bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() )
{
	return std::abs(x) < tolerance;
}