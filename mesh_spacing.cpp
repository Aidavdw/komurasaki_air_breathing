#include "mesh_spacing.h"
#include "field_quantity.h"
#include "domain.h"
#include <stdexcept>
#include <functional>
#include "gradient_descent.h"

MeshSpacing::MeshSpacing(const EMeshSpacingType meshSpacingType, const double elementLength, const int nElements, const double resolution_left, const double resolution_right) :
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
	bool bLeft = !IsCloseToZero(left);
	bool bRight = !IsCloseToZero(right);
	bool bN = (amountOfElements != 0);

	int valuesProvided = bLeft + bRight + bN;

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

	if (valuesProvided < (3 - requiredParameterCount))
		throw std::logic_error("The spacing type is underconstrained! More parameters need to be set.");
	if (valuesProvided > (3 - requiredParameterCount))
		throw std::logic_error("The spacing type is overconstrained! More parameters need to be left unset.");


	// Create a lambda with the given parameters always filled in, and the others as inputs
	std::vector<double> knowns = { left, right, float(amountOfElements) };
	auto glambda = [knowns](std::vector<double>* funcLoc)
	{
		int paramNumber = 0;
		std::vector<double> vec = knowns;
		for (int i = 0; i < knowns.size(); i++)
		{
			if (IsCloseToZero(knowns[i]))
			{
				vec[i] = funcLoc->at(paramNumber);
				paramNumber++;
			}
		}
		return SpacingObjectiveFunction(vec);
	};

	// Set up gradient descent
	qbGradient solver;

	// Assign the objective function.
	std::function<double(std::vector<double>*)> p_ObjectFcn{ glambda };
	solver.SetObjectFcn(p_ObjectFcn);

	// Set a start point.
	double knownSpacing = std::max({ left, right });
	double leftGuess = bLeft ? left : knownSpacing;
	double rightGuess = bRight ? right : knownSpacing;
	solver.SetStartPoint({ leftGuess, rightGuess, double(amountOfElements) });
	solver.SetMaxIterations(50);
	solver.SetStepSize(0.1);

	// Run gradient descent
	std::vector<double> funcLoc;
	double funcVal;
	solver.Optimize(&funcLoc, &funcVal);

}

bool IsCloseToZero(const double x, const double tolerance=std::numeric_limits<double>::epsilon() )
{
	return std::abs(x) < tolerance;
}

double SpacingObjectiveFunction(std::vector<double>& funcLoc)
{
	// Regardless of the actual function that is being optimised, there are 4 constraints, directly given by the input conditions
}
